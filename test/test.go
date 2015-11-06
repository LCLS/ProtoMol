package main

import (
	"bufio"
	"encoding/binary"
	"errors"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
)

var debug bool
var verbose bool
var parallel bool
var errorfailure bool

var vLogger *log.Logger

func main() {
	flag.BoolVar(&debug, "debug", false, "Enable debugging of ProtoMol execution")
	flag.BoolVar(&verbose, "verbose", false, "Enable verbose output")
	flag.BoolVar(&parallel, "parallel", false, "Enable parallel ProtoMol execution")
	flag.BoolVar(&errorfailure, "errorfailure", false, "Force failure of on a single test failure")
	single := flag.String("single", "", "Single test to run")
	flag.Parse()

	// Setup Verbose Logger
	if verbose {
		vLogger = log.New(os.Stdout, "", log.LstdFlags)
	} else {
		vLogger = log.New(ioutil.Discard, "", log.LstdFlags)
	}

	// Update Path Variable
	os.Setenv("PATH", ".:"+os.Getenv("PATH"))

	// Check if we can find ProtoMol
	if _, err := exec.LookPath("ProtoMol"); err != nil {
		log.Fatalln("ProtoMol executable could not be found")
	}

	if parallel {
		// Check if we can find mpirun
		if _, err := exec.LookPath("mpirun"); err != nil {
			log.Fatalln("MPI executable could not be found")
		}

		// Check if protomol was built with MPI
		out, _ := exec.Command("ProtoMol").CombinedOutput()
		if strings.Contains(string(out), "No MPI compilation") {
			log.Fatalln("Parallel tests unavailable. Please recompile ProtoMol with MPI support.")
		}
	}

	// Create Output Directory
	os.Mkdir("tests/output", 0755)

	// Find Tests
	tests, err := filepath.Glob("tests/*.conf")
	if err != nil {
		log.Fatalln(err)
	}

	if len(*single) > 0 {
		tests = []string{*single}
	}

	// Run Tests
	for _, test := range tests {
		RunTest(test, false)
		if parallel {
			RunTest(test, true)
		}
	}
}

func RunTest(config string, parallel bool) {
	basename := strings.Replace(filepath.Base(config), filepath.Ext(config), "", -1)

	var cmd *exec.Cmd
	if parallel {
		log.Println("Executing Parallel Test:", basename)
		cmd = exec.Command("mpirun", "-np", "2", "ProtoMol", config)
	} else {
		log.Println("Executing Test:", basename)
		cmd = exec.Command("ProtoMol", config)
	}

	if debug {
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stdout
	} else {
		cmd.Stdout = nil
		cmd.Stderr = nil
	}

	if err := cmd.Run(); err != nil {
		log.Println("\tIgnored")
		return
	}

	// Find Expected Results
	expects, err := filepath.Glob("tests/expected/" + basename + ".*")
	if err != nil {
		log.Fatalln(err)
	}

	// Compare Results
	for _, expected := range expects {
		output := "tests/output/" + filepath.Base(expected)
		extension := filepath.Ext(output)
		if extension == ".header" {
			continue
		}

		log.Println(fmt.Sprintf("\tTesting %s vs %s", output, expected))

		var result bool
		switch extension {
		case ".dcd":
			result = isMatchingDCD(output, expected)
			break
		case ".energy":
			result = isMatchingEnergy(output, expected)
			break
		case ".forces":
			result = isMatchingForce(output, expected)
			break
		case ".pos":
			result = isMatchingPosition(output, expected)
			break
		case ".vel":
			result = isMatchingVelocity(output, expected)
			break
		case ".vec":
			result = isMatchingVelocity(output, expected)
			break
		case ".val":
			result = isMatchingVelocity(output, expected)
			break
		default:
			log.Println("\t\tIgnored")
			continue
		}

		if !result {
			if errorfailure {
				log.Fatalln("\t\tFailed")
			} else {
				log.Println("\t\tFailed")
			}
		} else {
			log.Println("\t\tPassed")
		}
	}
}

// Generic Parsing
type AtomPosition struct {
	Name    string
	X, Y, Z float64
}

func ParseAtomPosition(line string) (*AtomPosition, error) {
	var values []string
	for _, value := range strings.Split(line, " ") {
		if len(value) > 0 {
			values = append(values, strings.TrimSpace(value))
		}
	}

	if len(values) != 4 {
		return nil, errors.New("invalid atom position line")
	}

	result := AtomPosition{Name: values[0]}

	x, err := strconv.ParseFloat(values[1], 64)
	if err != nil {
		return nil, err
	}
	result.X = x

	y, err := strconv.ParseFloat(values[2], 64)
	if err != nil {
		return nil, err
	}
	result.Y = y

	z, err := strconv.ParseFloat(values[3], 64)
	if err != nil {
		return nil, err
	}
	result.Z = z

	return &result, nil
}

// DCD Comparison
func isMatchingDCD(actual, expected string) bool {
	dcdActual, err := ReadDCD(actual)
	if err != nil {
		return true
	}

	dcdExpected, err := ReadDCD(expected)
	if err != nil {
		return true
	}

	if dcdActual.Header.Frames != dcdExpected.Header.Frames {
		vLogger.Println("Atom Count Differs")
		return false
	}

	if dcdActual.Header.FirstStep != dcdExpected.Header.FirstStep {
		vLogger.Printf("First Step Differs. Should be %d but is %d\n", dcdExpected.Header.FirstStep, dcdActual.Header.FirstStep)
		return false
	}

	if dcdActual.Atoms != dcdExpected.Atoms {
		vLogger.Printf("Atom Count Differs. Should be %d but is %d\n", dcdExpected.Atoms, dcdActual.Atoms)
		return false
	}

	diffs := 0
	for frame := int32(0); frame < dcdExpected.Header.Frames; frame++ {
		for atom := int32(0); atom < dcdExpected.Atoms; atom++ {
			ePosition := dcdExpected.Frames[frame].Position[atom]
			aPosition := dcdActual.Frames[frame].Position[atom]

			xDiff := math.Max(ePosition.X, aPosition.X) - math.Min(ePosition.X, aPosition.X)
			if xDiff > 0.00001 {
				diffs++
				vLogger.Printf("Frame %d, Atom %d Differs X. Expected: %f, Actual: %f, Difference: %f\n", frame, atom, ePosition.X, aPosition.X, xDiff)
			}

			yDiff := math.Max(ePosition.Y, aPosition.Y) - math.Min(ePosition.Y, aPosition.Y)
			if yDiff > 0.00001 {
				diffs++
				vLogger.Printf("Frame %d, Atom %d Differs Y. Expected: %f, Actual: %f, Difference: %f\n", frame, atom, ePosition.Y, aPosition.Y, yDiff)
			}

			zDiff := math.Max(ePosition.X, aPosition.X) - math.Min(ePosition.X, aPosition.X)
			if zDiff > 0.00001 {
				diffs++
				vLogger.Printf("Frame %d, Atom %d Differs Z. Expected: %f, Actual: %f, Difference: %f\n", frame, atom, ePosition.Z, aPosition.Z, zDiff)
			}
		}
	}

	return diffs == 0
}

// DCD Parsing
type Vector3f struct {
	X, Y, Z float64
}
type DCDFile struct {
	Header  DCDHeader
	Comment []byte
	Atoms   int32
	Frames  []DCDFrame
}

type DCDFrame struct {
	Position []Vector3f
}
type DCDHeader struct {
	Frames      int32
	FirstStep   int32
	_           [24]byte
	FreeIndexes int32
	_           [44]byte
	_           [4]byte
}

func ReadDCD(filename string) (*DCDFile, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	// Read Endian
	var endian int32
	if err := binary.Read(file, binary.LittleEndian, &endian); err != nil {
		return nil, err
	}

	// Test Endian
	var order binary.ByteOrder = binary.LittleEndian
	if endian != 84 {
		order = binary.BigEndian
	}

	// Read Cord
	var cord [4]byte
	if err := binary.Read(file, order, &cord); err != nil {
		return nil, err
	}

	// Check Magic Number
	if cord != [4]byte{'C', 'O', 'R', 'D'} {
		return nil, err
	}

	var result DCDFile

	// Read Header
	if err := binary.Read(file, order, &result.Header); err != nil {
		return nil, err
	}

	// Read Comment
	if err := FortranReadString(file, order, &result.Comment); err != nil {
		return nil, err
	}

	// Read Atoms
	if err := FortranRead(file, order, &result.Atoms); err != nil {
		return nil, err
	}

	// Skip Free Indexies
	if result.Header.FreeIndexes > 0 {
		io.CopyN(ioutil.Discard, file, int64(4*(result.Atoms-result.Header.FreeIndexes+2)))
	}

	// Read Frames
	tempX := make([]float32, result.Atoms)
	tempY := make([]float32, result.Atoms)
	tempZ := make([]float32, result.Atoms)

	for i := int32(0); i < result.Header.Frames; i++ {
		var frame DCDFrame

		// Read X
		if err := FortranRead(file, order, &tempX); err != nil {
			return nil, err
		}

		// Read Y
		if err := FortranRead(file, order, &tempY); err != nil {
			return nil, err
		}

		// Read Z
		if err := FortranRead(file, order, &tempZ); err != nil {
			return nil, err
		}

		for j := int32(0); j < result.Atoms; j++ {
			frame.Position = append(frame.Position, Vector3f{X: float64(tempX[j]), Y: float64(tempY[j]), Z: float64(tempZ[j])})
		}
		result.Frames = append(result.Frames, frame)
	}

	return &result, nil
}

func FortranRead(file io.Reader, order binary.ByteOrder, data interface{}) error {
	var head int32
	if err := binary.Read(file, order, &head); err != nil {
		return err
	}

	if err := binary.Read(file, order, data); err != nil {
		return err
	}

	var tail int32
	if err := binary.Read(file, order, &tail); err != nil {
		return err
	}

	if head != tail {
		return errors.New(fmt.Sprintf("head (%d) does not match tail (%d)", head, tail))
	}

	return nil
}

func FortranReadString(file io.Reader, order binary.ByteOrder, data *[]byte) error {
	var head int32
	if err := binary.Read(file, order, &head); err != nil {
		return err
	}

	*data = make([]byte, head)
	if err := binary.Read(file, order, data); err != nil {
		return err
	}

	var tail int32
	if err := binary.Read(file, order, &tail); err != nil {
		return err
	}

	return nil
}

// Energy Comparison
func isMatchingEnergy(actual, expected string) bool {
	energyActual, err := ReadEnergy(actual)
	if err != nil {
		return true
	}

	energyExpected, err := ReadEnergy(expected)
	if err != nil {
		return true
	}

	if len(energyActual.Column) != len(energyExpected.Column) {
		vLogger.Printf("Column Count Differs. Should be %d but is %d\n", len(energyActual.Column), len(energyExpected.Column))
		return false
	}

	var diffs int
	for i := 0; i < len(energyExpected.Column); i++ {
		if energyActual.Column[i].Name != energyExpected.Column[i].Name {
			vLogger.Printf("Column Name Differs. Should be %s but is %s\n", energyActual.Column[i].Name, energyExpected.Column[i].Name)
			return false
		}

		if len(energyActual.Column[i].Values) != len(energyExpected.Column[i].Values) {
			vLogger.Printf("Column Value Count Differs. Should be %d but is %d\n", len(energyActual.Column[i].Values), len(energyExpected.Column[i].Values))
			return false
		}

		for j := 0; j < len(energyExpected.Column[i].Values); j++ {
			xExpected := energyExpected.Column[i].Values[j]
			xActual := energyActual.Column[i].Values[j]

			if math.Max(float64(xExpected), float64(xActual))-math.Min(float64(xExpected), float64(xActual)) > 0.00001 {
				diffs++
				vLogger.Printf("Column Value Differs. Expected: %f, Actual: %f, Difference: %f\n", xExpected, xActual, math.Max(float64(xExpected), float64(xActual))-math.Min(float64(xExpected), float64(xActual)))
			}
		}
	}

	return diffs == 0
}

// Energy Parsing
type EnergyFile struct {
	Column []EnergyFileColumn
}

type EnergyFileColumn struct {
	Name   string
	Values []float32
}

func ReadEnergy(filename string) (*EnergyFile, error) {
	// Parse Headers
	sHeader, err := ioutil.ReadFile(fmt.Sprintf("%s.header", filename))
	if err != nil {
		return nil, err
	}

	var result EnergyFile

	for _, header := range strings.Split(string(sHeader), " ") {
		if len(header) == 0 {
			continue
		}
		result.Column = append(result.Column, EnergyFileColumn{Name: header})
	}

	// Parse File
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		columnID := 0
		for _, column := range strings.Split(scanner.Text(), " ") {
			if len(column) == 0 {
				continue
			}

			value, err := strconv.ParseFloat(column, 32)
			if err != nil {
				return nil, err
			}

			result.Column[columnID].Values = append(result.Column[columnID].Values, float32(value))
			columnID++
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return &result, nil
}

// Force Comparison
func isMatchingForce(actual, expected string) bool {
	aForce, err := ReadForce(actual)
	if err != nil {
		return true
	}

	eForce, err := ReadForce(expected)
	if err != nil {
		return true
	}

	if aForce.FrameCount != eForce.FrameCount {
		vLogger.Printf("Frame Count Differs. Should be %d but is %d\n", eForce.FrameCount, aForce.FrameCount)
		return false
	}

	if aForce.AtomCount != eForce.AtomCount {
		vLogger.Printf("Atom Count Differs. Should be %d but is %d\n", eForce.AtomCount, aForce.AtomCount)
		return false
	}

	diffs := 0
	for frame := 0; frame < eForce.FrameCount; frame++ {
		for atom := 0; atom < eForce.AtomCount; atom++ {
			xExpected := eForce.Frame[frame].Atom[atom].X
			xActual := aForce.Frame[frame].Atom[atom].X

			if math.Max(xExpected, xActual)-math.Min(xExpected, xActual) > 0.00001 {
				diffs++
				vLogger.Printf("Frame %d, Atom %d Differs. Expected: %f, Actual: %f, Difference: %f\n", frame, atom, xExpected, xActual, math.Max(xExpected, xActual)-math.Min(xExpected, xActual))
			}

			yExpected := eForce.Frame[frame].Atom[atom].Y
			yActual := aForce.Frame[frame].Atom[atom].Y

			if math.Max(yExpected, yActual)-math.Min(yExpected, yActual) > 0.00001 {
				diffs++
				vLogger.Printf("Frame %d, Atom %d Differs. Expected: %f, Actual: %f, Difference: %f\n", frame, atom, yExpected, yActual, math.Max(yExpected, yActual)-math.Min(yExpected, yActual))
			}

			zExpected := eForce.Frame[frame].Atom[atom].Z
			zActual := aForce.Frame[frame].Atom[atom].Z

			if math.Max(zExpected, zActual)-math.Min(zExpected, zActual) > 0.00001 {
				diffs++
				vLogger.Printf("Frame %d, Atom %d Differs. Expected: %f, Actual: %f, Difference: %f\n", frame, atom, zExpected, zActual, math.Max(zExpected, zActual)-math.Min(zExpected, zActual))
			}
		}
	}

	return diffs == 0
}

// Energy Parsing
type ForceFile struct {
	FrameCount int
	AtomCount  int
	Frame      []ForceFileFrame
}

type ForceFileFrame struct {
	Atom []AtomPosition
}

func ReadForce(filename string) (*ForceFile, error) {
	// Parse File
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	r := bufio.NewReader(file)

	// Read Framecount and atomcount
	sFrames, err := r.ReadString('\n')
	if err != nil {
		return nil, err
	}

	sAtoms, err := r.ReadString('\n')
	if err != nil {
		return nil, err
	}

	// Convert framecount and atomcount
	frameCount, err := strconv.Atoi(sFrames)
	if err != nil {
		return nil, err
	}

	atomCount, err := strconv.Atoi(sAtoms)
	if err != nil {
		return nil, err
	}

	result := ForceFile{FrameCount: frameCount, AtomCount: atomCount}

	// Parse Frames
	for frame := 0; frame < result.FrameCount; frame++ {
		var frameData ForceFileFrame
		for atom := 0; atom < result.AtomCount; atom++ {
			sAtom, err := r.ReadString('\n')
			if err != nil {
				return nil, err
			}

			atomValue, err := ParseAtomPosition(sAtom)
			if err != nil {
				return nil, err
			}

			frameData.Atom = append(frameData.Atom, *atomValue)
		}
		result.Frame = append(result.Frame, frameData)
	}

	return &result, nil
}

// Position Comparison
func isMatchingPosition(actual, expected string) bool {
	aPosition, err := ReadPosition(actual)
	if err != nil {
		return true
	}

	ePosition, err := ReadPosition(expected)
	if err != nil {
		return true
	}

	if aPosition.AtomCount != ePosition.AtomCount {
		vLogger.Printf("Atom Count Differs. Should be %d but is %d\n", ePosition.AtomCount, aPosition.AtomCount)
		return false
	}

	diffs := 0
	for atom := 0; atom < ePosition.AtomCount; atom++ {
		xExpected := ePosition.Atom[atom].X
		xActual := aPosition.Atom[atom].X

		if math.Max(xExpected, xActual)-math.Min(xExpected, xActual) > 0.00001 {
			diffs++
			vLogger.Printf("Atom %d Differs. Expected: %f, Actual: %f, Difference: %f\n", atom, xExpected, xActual, math.Max(xExpected, xActual)-math.Min(xExpected, xActual))
		}

		yExpected := ePosition.Atom[atom].Y
		yActual := aPosition.Atom[atom].Y

		if math.Max(yExpected, yActual)-math.Min(yExpected, yActual) > 0.00001 {
			diffs++
			vLogger.Printf("Atom %d Differs. Expected: %f, Actual: %f, Difference: %f\n", atom, yExpected, yActual, math.Max(yExpected, yActual)-math.Min(yExpected, yActual))
		}

		zExpected := ePosition.Atom[atom].Z
		zActual := aPosition.Atom[atom].Z

		if math.Max(zExpected, zActual)-math.Min(zExpected, zActual) > 0.00001 {
			diffs++
			vLogger.Printf("Atom %d Differs. Expected: %f, Actual: %f, Difference: %f\n", atom, zExpected, zActual, math.Max(zExpected, zActual)-math.Min(zExpected, zActual))
		}
	}

	return diffs == 0
}

// Energy Parsing
type PositionFile struct {
	AtomCount int
	Atom      []AtomPosition
}

func ReadPosition(filename string) (*PositionFile, error) {
	// Parse File
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	r := bufio.NewReader(file)

	// Read AtomCount and convert
	sAtoms, err := r.ReadString('\n')
	if err != nil {
		return nil, err
	}

	atomCount, err := strconv.Atoi(sAtoms)
	if err != nil {
		return nil, err
	}

	result := PositionFile{AtomCount: atomCount}

	// Parse Frames
	for atom := 0; atom < result.AtomCount; atom++ {
		sAtom, err := r.ReadString('\n')
		if err != nil {
			return nil, err
		}

		atomValue, err := ParseAtomPosition(sAtom)
		if err != nil {
			return nil, err
		}

		result.Atom = append(result.Atom, *atomValue)
	}

	return &result, nil
}

// Velocity Comparison
func isMatchingVelocity(actual, expected string) bool {
	aPosition, err := ReadPosition(actual)
	if err != nil {
		return true
	}

	ePosition, err := ReadPosition(expected)
	if err != nil {
		return true
	}

	if aPosition.AtomCount != ePosition.AtomCount {
		vLogger.Printf("Atom Count Differs. Should be %d but is %d\n", ePosition.AtomCount, aPosition.AtomCount)
		return false
	}

	diffs := 0
	for atom := 0; atom < ePosition.AtomCount; atom++ {
		xExpected := ePosition.Atom[atom].X
		xActual := aPosition.Atom[atom].X

		if math.Max(xExpected, xActual)-math.Min(xExpected, xActual) > 0.00001 {
			diffs++
			vLogger.Printf("Atom %d Differs. Expected: %f, Actual: %f, Difference: %f\n", atom, xExpected, xActual, math.Max(xExpected, xActual)-math.Min(xExpected, xActual))
		}

		yExpected := ePosition.Atom[atom].Y
		yActual := aPosition.Atom[atom].Y

		if math.Max(yExpected, yActual)-math.Min(yExpected, yActual) > 0.00001 {
			diffs++
			vLogger.Printf("Atom %d Differs. Expected: %f, Actual: %f, Difference: %f\n", atom, yExpected, yActual, math.Max(yExpected, yActual)-math.Min(yExpected, yActual))
		}

		zExpected := ePosition.Atom[atom].Z
		zActual := aPosition.Atom[atom].Z

		if math.Max(zExpected, zActual)-math.Min(zExpected, zActual) > 0.00001 {
			diffs++
			vLogger.Printf("Atom %d Differs. Expected: %f, Actual: %f, Difference: %f\n", atom, zExpected, zActual, math.Max(zExpected, zActual)-math.Min(zExpected, zActual))
		}
	}

	return diffs == 0
}
