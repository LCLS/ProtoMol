#!/usr/bin/python

import sys
import os
import glob
import shlex
import subprocess
import comparator
import compare_dcd
import argparse

import logging

DEFAULT_EPSILON = 0.00001
DEFAULT_SCALINGFACTOR = 1.0

def parse_params(flname):
    """
....Parses test parameters from Protomol configuration files.
....Paramaters are used to change how testing is performed.
....It is assumed that parameters will be in the form:
            ## key = value
....where key will be a string, value will by any python structure
....that fits on one line, and the hashmarks will be the first two
        characters on the line.
...."""

    fl = open(flname)
    params = dict()
    for ln in fl:

        if ln[:2] != '##':
            continue
        eq_pos = ln.find('=')
        if eq_pos == -1:
            continue
        param_name = ln[2:eq_pos].strip()
        param_str = ln[eq_pos + 1:].strip()
        param_val = eval(param_str)
        params[param_name] = param_val
    fl.close()
    return params


def run_test(protomol_path, conf_file, pwd, parallel):
    tests = 0
    testspassed = 0
    testsfailed = 0
    failedtests = []

    conf_param_overrides = parse_params(conf_file)
    epsilon = conf_param_overrides.get('epsilon', DEFAULT_EPSILON)
    scaling_factor = conf_param_overrides.get('scaling_factor', DEFAULT_SCALINGFACTOR)

    base = os.path.splitext(os.path.basename(conf_file))[0]
    if not parallel:
        logging.info('Executing Test: ' + base)
    else:
        logging.info('Executing Parallel Test: ' + base)

    cmd = []
    if parallel:
        cmd.append('mpirun')
        cmd.append('-np')
        cmd.append('2')
    cmd.append(protomol_path)
    cmd.append(conf_file)

    DEVNULL = open(os.devnull, 'w')
    p = subprocess.Popen(cmd, stdout=None, stderr=None)
    p.communicate()
    if p.returncode > 0:
        s = 'Not able to execute Protomol!\n'
        s += 'cmd: ' + str(cmd) + '\n'
        logging.critical(s)
    expects = []
    outputs = glob.glob('tests/output/' + base + '.*')
    outputtemp = []
    for output in outputs:
        outbase = os.path.basename(output)
        if os.path.exists('tests/expected/' + outbase):
            outputtemp.append(output)
            expects.append('tests/expected/' + outbase)

    outputs = outputtemp

    for i in xrange(len(outputs)):
        tests += 1
        ftype = os.path.splitext(os.path.basename(outputs[i]))[1]

        if ftype in ['.header', '.xtc']:
            continue

        ignoreSign = False

        # Ignore signs on eignevectors
        if ftype == '.vec':
            ignoreSign = True

        logging.info('\tTesting: ' + expects[i] + ' ' + outputs[i])

        if ftype == ".dcd":
            if compare_dcd.compare_dcd(expects[i], outputs[i], epsilon, scaling_factor, ignoreSign):
                logging.info('\t\tPassed')
                testspassed += 1
            else:
                logging.warning('\t\tFailed')
                testsfailed += 1
                failedtests.append('Comparison of ' + expects[i] + ' and ' + outputs[i])
                if args.errorfailure:
                    sys.exit(1)
        else:
            if comparator.compare(expects[i], outputs[i], epsilon, scaling_factor, ignoreSign):
                logging.info('\t\tPassed')
                testspassed += 1
            else:
                logging.warning('\t\tFailed')
                testsfailed += 1
                failedtests.append('Comparison of ' + expects[i] + ' and ' + outputs[i])
                if args.errorfailure:
                    sys.exit(1)

    return (tests, testspassed, testsfailed, failedtests)


def find_conf_files(pwd, args):
    files = []
    if args.single == None and args.regex == None:
        files = glob.glob(os.path.join(pwd, 'tests', '*.conf'))
    elif args.single != None:
        if args.single.find('.conf') != -1:
            files.append(args.single)
        else:
            raise Exception, args.single + ' is not a valid configuration file'
    elif args.regex != None:
        files = glob.glob(os.path.join(pwd, 'tests', args.regex + '.conf'))
    return files


def find_protomol(pwd):
    if os.path.exists(os.path.join(pwd, 'ProtoMol')):
        return os.path.join(pwd, 'ProtoMol')
    elif os.path.exists(os.path.join(pwd, 'ProtoMol.exe')):
        return os.path.join(pwd, 'ProtoMol.exe')
    else:
        raise Exception, 'Cannot find ProtoMol executable in ' + pwd + ' . Please put the ProtoMol executable in this directory'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ProtoMol Test Suite')
    parser.add_argument('--verbose', '-v', action='store_true', default=False, help='Verbose output')
    parser.add_argument('--errorfailure', '-e', action='store_true', default=False, help='Break on test failure')
    parser.add_argument('--serial', '-s', action='store_true', default=True, help='Serial Testing')
    parser.add_argument('--parallel', '-p', action='store_true', default=False, help='MPI Testing')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--single', '-t', help='Single test to run. Must be within the tests directory.')
    group.add_argument('--regex', '-r', help='Regular expression of tests to run, Requires quotation marks around argument')

    args = parser.parse_args()

    level = logging.INFO
    if args.verbose:
        level = logging.DEBUG
    logging.basicConfig(level=level, format='%(message)s')

    pwd = os.getcwd()
    files = find_conf_files(pwd, args)
    files.sort()
    logging.debug('Files: ' + str(files))

    protomol_path = find_protomol(pwd)
    if args.parallel:
        p = subprocess.Popen(protomol_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        if "No MPI compilation" in stderr:
            logging.critical("Parallel tests unavailable. Please recompile ProtoMol with MPI support.")
            sys.exit(1)

    tests = 0
    testspassed = 0
    testsfailed = 0
    failedtests = []

    for conf_file in files:
        if args.serial:
            results = run_test(protomol_path, conf_file, pwd, False)
            tests += results[0]
            testspassed += results[1]
            testsfailed += results[2]
            failedtests.extend(results[3])
        if args.parallel:
            results = run_test(protomol_path, conf_file, pwd, True)
            tests += results[0]
            testspassed += results[1]
            testsfailed += results[2]
            failedtests.extend(results[3])

    testsnotrun = tests - (testspassed + testsfailed)

    logging.info('')
    logging.info('Tests: %d' % tests)
    logging.info('Tests Not Run: %d' % testsnotrun)
    logging.info('')
    logging.info('Tests Passed: %d' % testspassed)
    logging.info('Tests Failed: %d\n' % testsfailed)

    if len(failedtests) > 0:
        for test in failedtests:
            logging.warning('%s failed.' % test)
        sys.exit(1)
