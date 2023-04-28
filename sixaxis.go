// Package sixaxis does forward and inverse kinematics for a six
// axis robot defined as a sequence of radial joints.
//
// The radial joints are arranged in the following global coordinates
// with the robot in its "all angles zero" pose (all joints point arm
// segment up in positive Z direction).
//
//	J0 = rotate around Z axis an arm segment of length d0
//	J1 = rotate around X axis an arm segment of length d1
//	J2 = rotate around X axis an arm segment of length d2
//	J3 = rotate around Z axis an arm segment of length d3
//	J4 = rotate around X axis an arm segment of length d4
//	J5 = rotate around Z axis an arm segment of length d5 = 0
package sixaxis

import (
	"errors"
	"fmt"
	"math"

	"zappem.net/pub/math/geom"
)

// Joint holds a robot joint configuration.
type Joint struct {
	Min, Max geom.Angle
	Rot      func(a geom.Angle) geom.Matrix
	Offset   float64
}

// Pose holds the current pose of the robot in both cartesian and
// joint coordinates.
type Pose struct {
	M geom.Matrix
	V geom.Vector
	J []geom.Angle
}

// Robot holds the combined geometric state of a robot.
type Robot struct {
	// Precision holds the dtheta precision of the robot. It
	// defaults to geom.Zeroish fraction of a full rotation.
	Precision geom.Angle

	j []Joint
	p Pose
}

// Pose returns the current pose of the robot.
func (r *Robot) Pose() Pose {
	return r.p
}

// NewRobot specifies the joints for a robot. Its default pose is all
// joint parameters are zero.
func NewRobot(js ...Joint) (*Robot, error) {
	if len(js) != 6 {
		return nil, fmt.Errorf("require six joint parameters not %d", len(js))
	}
	r := &Robot{
		j: js,
		p: Pose{
			J: make([]geom.Angle, len(js)),
		},
		Precision: geom.Degrees(360 * geom.RefZero()),
	}
	r.fwd()
	return r, nil
}

// J returns the joint parameter for the specified joint.
func (r *Robot) J(i int) geom.Angle {
	if i < 0 || i > len(r.j) {
		return 0
	}
	return r.p.J[i]
}

// MV returns the orientation and offset of the J5 face-plate of
// robot. This is equal to the last forward kinematics solution for
// the robot's pose.
func (r *Robot) MV() (geom.Matrix, geom.Vector) {
	return r.p.M, r.p.V
}

// Forward evaluates the forward kinematics for a set of joint angles and
// returns the basis and vector of the face-plate. It does not alter the
// working pose of the robot.
func (r *Robot) Forward(j []geom.Angle) (geom.Matrix, geom.Vector) {
	m := geom.I
	v := geom.V(0)
	for i := len(r.j) - 1; i >= 0; i-- {
		joint := r.j[i].Rot(j[i])
		v = joint.XV(v.Add(geom.Z(r.j[i].Offset)))
		m = joint.XM(m)
	}
	return m, v
}

// fwd computes the forward kinematics given the currently known joint
// angles and sets these as the current pose of the robot. Use this
// function to fixup the robot's pose after changing its joints.
func (r *Robot) fwd() error {
	r.p.M, r.p.V = r.Forward(r.p.J)
	return nil
}

// Err* are the errors exported by this package.
var (
	ErrBadJoint   = errors.New("invalid joint")
	ErrLimit      = errors.New("paramter outside joint range")
	ErrNoSolution = errors.New("no inverse kinimatics solution")
)

// SetJ forces the joint of the robot to a joint parameter value. An
// error is returned if the joint cannot adopt that paramter because
// of a range limit.
func (r *Robot) SetJ(i int, a geom.Angle) error {
	if i < 0 || i > len(r.j) {
		return ErrBadJoint
	}
	if j := r.j[i]; a <= j.Min || a >= j.Max {
		return ErrLimit
	}
	r.p.J[i] = a
	return r.fwd()
}

// validate compares a set of angles with the supplied pose of the
// robot and should the angles generate the pose, they are forwarded
// to soln.
func (r *Robot) validate(soln chan<- []geom.Angle, b geom.Matrix, v geom.Vector, a ...geom.Angle) {
	j := r.j[5]
	if j.Min != j.Max && (a[5] < j.Min || a[5] > j.Max) {
		// This angle is out of range so no solution here.
		return
	}
	m, u := r.Forward(a)
	if !m.Equals(b) || !u.Equals(v) {
		return // no match
	}
	soln <- a

	// Allow for 360 deg rotations too - j5 could be like that.
	a5 := a[5] + 2*math.Pi
	if a5 >= 0 {
		a5 = a[5] - 2*math.Pi
		if a5 <= j.Min {
			return
		}
	} else if a5 >= j.Max {
		return
	}

	soln <- []geom.Angle{a[0], a[1], a[2], a[3], a[4], a5}
}

func (r *Robot) solve5(soln chan<- []geom.Angle, b geom.Matrix, v geom.Vector, a0, a1, a2, a3, a4 geom.Angle) {
	if j := r.j[3]; j.Min != j.Max && (a3 <= j.Min || a3 >= j.Max) {
		// This angle is out of range so no solution here.
		return
	}

	B00 := b[0]
	B10 := b[3]
	B20 := b[6]
	s0 := a0.S()
	c0 := a0.C()
	a12 := a1 + a2
	s12 := a12.S()
	c12 := a12.C()
	s4 := a4.S()

	if !geom.Zeroish(s4) {
		s5 := ((B00*s0-B10*c0)*s12 + B20*c12) / s4
		if s5*s5 > 1 {
			// No solution.
			return
		}
		j5 := geom.Angle(math.Asin(s5))

		r.validate(soln, b, v, a0, a1, a2, a3, a4, j5)
		if j5p := j5.LikeSin(); !geom.Zeroish(float64(j5 - j5p)) {
			r.validate(soln, b, v, a0, a1, a2, a3, a4, j5p)
		}
		return
	}

	s3 := a3.S()
	B01 := b[1]
	B11 := b[4]

	x := -B00*c0 - B10*s0
	y := -B01*c0 - B11*s0

	if geom.Zeroish(s3) {
		j5 := r.p.J[5]
		// preserve current j5 if unable to limit value.
		if !geom.Zeroish(x) || !geom.Zeroish(y) {
			j5 = geom.Angle(math.Atan2(-y, x))
			if j5p := j5.LikeTan(); !geom.Zeroish(float64(j5 - j5p)) {
				r.validate(soln, b, v, a0, a1, a2, a3, a4, j5p)
			}
		}
		r.validate(soln, b, v, a0, a1, a2, a3, a4, j5)
		return
	}

	c3 := a3.C()
	c4 := a4.C()

	s5 := (x*c4*s3 + c3*y) / (c4*s3 + c3*c3)
	if s5*s5 > 1 {
		return
	}

	j5 := geom.Angle(math.Asin(s5))
	r.validate(soln, b, v, a0, a1, a2, a3, a4, j5)
	if j5p := j5.LikeSin(); !geom.Zeroish(float64(j5 - j5p)) {
		r.validate(soln, b, v, a0, a1, a2, a3, a4, j5p)
	}
}

func (r *Robot) solve3(soln chan<- []geom.Angle, b geom.Matrix, v geom.Vector, a0, a1, a2, a4 geom.Angle) {
	if j := r.j[4]; j.Min != j.Max && (a4 <= j.Min || a4 >= j.Max) {
		// This angle is out of range so no solution here.
		return
	}

	d4 := r.j[4].Offset
	r0 := v[0]
	r1 := v[1]
	s4 := a4.S()

	if !geom.Zeroish(s4) {
		s3 := (a0.C()*r0 + a0.S()*r1) / (s4 * d4)
		if s3*s3 > 1 {
			// No solution.
			return
		}
		j3 := geom.Angle(math.Asin(s3))
		r.solve5(soln, b, v, a0, a1, a2, j3, a4)
		if j3p := j3.LikeSin(); !geom.Zeroish(float64(j3 - j3p)) {
			r.solve5(soln, b, v, a0, a1, a2, j3p, a4)
		}
		return
	}

	B02 := b[2]
	B12 := b[5]
	B22 := b[8]
	c0 := a0.C()
	s0 := a0.S()
	a12 := a1 + a2
	c12 := a12.C()
	s12 := a12.S()

	x := B02*c0 + B12*s0
	y := (B12*c0-B02*s0)*c12 + B22*s12
	j3 := r.p.J[3]
	if !geom.Zeroish(x) || !geom.Zeroish(y) {
		j3 = geom.Angle(math.Atan2(y, x))
		if j3p := j3.LikeTan(); !geom.Zeroish(float64(j3 - j3p)) {
			r.solve5(soln, b, v, a0, a1, a2, j3p, a4)
		}
	}
	r.solve5(soln, b, v, a0, a1, a2, j3, a4)
}

func (r *Robot) solve4(soln chan<- []geom.Angle, b geom.Matrix, v geom.Vector, a0, a1, a2 geom.Angle) {
	if j := r.j[2]; j.Min != j.Max && (a2 <= j.Min || a2 >= j.Max) {
		// This angle is out of range so no solution here.
		return
	}

	a12 := a1 + a2
	s12 := a12.S()
	c12 := a12.C()
	B02 := b[2]
	B12 := b[5]
	B22 := b[8]

	c4 := s12*(B02*a0.S()-B12*a0.C()) + B22*c12
	if c4*c4 > 1 {
		// No solution.
		return
	}
	j4 := geom.Angle(math.Acos(c4))
	r.solve3(soln, b, v, a0, a1, a2, j4)
	if j4p := j4.LikeCos(); !geom.Zeroish(float64(j4 - j4p)) {
		r.solve3(soln, b, v, a0, a1, a2, j4p)
	}
}

func (r *Robot) solve2(soln chan<- []geom.Angle, b geom.Matrix, v geom.Vector, a0, a1 geom.Angle, x, y float64) {
	if j := r.j[1]; j.Min != j.Max && (a1 <= j.Min || a1 >= j.Max) {
		// This angle is out of range so no solution here.
		return
	}
	d1 := r.j[1].Offset
	d2 := r.j[2].Offset
	d3 := r.j[3].Offset

	s12 := (x - d1*a1.S()) / (d2 + d3)
	c12 := (y - d1*a1.C()) / (d2 + d3)

	if s12*s12 > 1 || c12*c12 > 1 {
		// No valid solution.
		return
	}

	// default to not changing the j2 angle if it can't be constrained.
	j2 := r.p.J[2]
	if !geom.Zeroish(s12) || !geom.Zeroish(c12) {
		a12 := geom.Angle(math.Atan2(s12, c12))
		j2 = a12 - a1
	}
	r.solve4(soln, b, v, a0, a1, j2)
}

func (r *Robot) solve1(soln chan<- []geom.Angle, b geom.Matrix, v geom.Vector, a0 geom.Angle) {
	if j := r.j[0]; j.Min != j.Max && (a0 <= j.Min || a0 >= j.Max) {
		// This angle is out of range so no solution here.
		return
	}
	r0 := v[0]
	r1 := v[1]
	r2 := v[2]
	d0 := r.j[0].Offset
	d1 := r.j[1].Offset
	d2 := r.j[2].Offset
	d3 := r.j[3].Offset
	d4 := r.j[4].Offset
	d23 := d2 + d3
	b02 := b[2]
	b12 := b[5]
	b22 := b[8]
	s0 := a0.S()
	c0 := a0.C()

	x := (-b02*s0+b12*c0)*d4 - c0*r1 + r0*s0
	y := -b22*d4 - d0 + r2
	x2y2 := x*x + y*y
	q := x2y2 + d1*d1 - d23*d23

	if geom.Zeroish(x2y2) {
		return // no solution.
	}
	if geom.Zeroish(x) {
		c1 := q / (2 * d1 * y)
		if c1*c1 > 1 {
			return // no solution
		}
		j1 := geom.Angle(math.Acos(c1))
		r.solve2(soln, b, v, a0, j1, x, y)
		if j1p := j1.LikeCos(); !geom.Zeroish(float64(j1 - j1p)) {
			r.solve2(soln, b, v, a0, j1p, x, y)
		}
		return

	}
	if geom.Zeroish(y) {
		s1 := q / (2 * d1 * x)
		if s1*s1 > 1 {
			return // no solution
		}
		j1 := geom.Angle(math.Asin(s1))
		r.solve2(soln, b, v, a0, j1, x, y)
		r.solve2(soln, b, v, a0, j1.LikeSin(), x, y)
		return
	}

	sq := 4*d1*d1*x2y2 - q*q
	if geom.Zeroish(sq) {
		sq = 0
	} else if sq < 0 {
		return
	}

	rt := math.Sqrt(sq)
	s1 := (q*x + y*rt) / (2 * d1 * x2y2)
	if s1*s1 <= 1 {
		j1 := geom.Angle(math.Asin(s1))
		r.solve2(soln, b, v, a0, j1, x, y)
		if j1p := j1.LikeSin(); !geom.Zeroish(float64(j1 - j1p)) {
			r.solve2(soln, b, v, a0, j1p, x, y)
		}
	}
	if geom.Zeroish(sq) {
		return
	}
	s1 = (q*x - y*rt) / (2 * d1 * x2y2)
	if s1*s1 <= 1 {
		j1 := geom.Angle(math.Asin(s1))
		r.solve2(soln, b, v, a0, j1, x, y)
		if j1p := j1.LikeSin(); !geom.Zeroish(float64(j1 - j1p)) {
			r.solve2(soln, b, v, a0, j1p, x, y)
		}
	}
}

// ikSolve computes all of the solutions to the inverse kinematics for
// the robot and returns them all via the soln channel.
func (r *Robot) ikSolve(soln chan<- []geom.Angle, b geom.Matrix, v geom.Vector) {
	defer close(soln)
	r0 := v[0]
	r1 := v[1]
	d4 := r.j[4].Offset
	b02 := b[2]
	b12 := b[5]
	x := r1 - b12*d4
	y := b02*d4 - r0

	// default to not changing the j0 angle if it can't be constrained.
	j0 := r.p.J[0]
	if geom.Zeroish(x) && geom.Zeroish(y) {
		r.solve1(soln, b, v, j0)
		return
	}

	j0 = geom.Angle(math.Atan2(y, x))
	r.solve1(soln, b, v, j0)
	r.solve1(soln, b, v, j0.LikeTan())
}

// Inverse computes the potential set of inverse kinematics for a
// target RV pose. It is the responsibility of the client to pick
// which of the solutions to use. If no angle selections are returned
// there are no solutions.
func (r *Robot) Inverse(m geom.Matrix, v geom.Vector) [][]geom.Angle {
	soln := make(chan []geom.Angle)
	go r.ikSolve(soln, m, v)
	var ans [][]geom.Angle
	for s := range soln {
		ans = append(ans, s)
	}
	return ans
}

// Closest picks the closest pose to the ref (or if ref is nil the
// current) pose from a list of equivalent joint positions. This
// function currently treats all joints with equal weight.
func (r *Robot) Closest(ref *Pose, cjs [][]geom.Angle) (int, error) {
	if ref == nil {
		ref = &r.p
	}
	best := -1
	var bestA2 geom.Angle
	for i, js := range cjs {
		var a2 geom.Angle
		for j, a := range js {
			d := a - ref.J[j]
			a2 += d * d
		}
		if best == -1 || bestA2 > a2 {
			best = i
			bestA2 = a2
		}
	}
	if best != -1 {
		return best, nil
	}
	return best, ErrNoSolution
}

// Pace represents a fraction completed of a linear path between two
// poses in terms of a number [0,1] and a collection of joint angles
// for the robot. Typically, Paces are generated in a way that
// interpolating all joints between two paces equally will follow the
// desired trajectory.
type Pace struct {
	Frac float64
	J    []geom.Angle
}

// linear uses binary search to find a series of succesive waypoints
// (expressed as Paces) between was and target (as parameterized by an
// angle, b in (0,a]) that satisfies the constraint that all joints
// appear to be linearly interpolated between successive waypoints to
// b. The function returns this set of joints and the fraction of the
// path completed, at the completion of each Pace, over the Paces
// channel. Linearity of interpolation is determined by computing the
// mid point and validating that its joint angles are 50% between the
// end point joint angles. The accuracy of the 50% test is in terms of
// the robot's Precision number. That is, the cumulative error from
// the rotations on each joint amounts to a displacement in x-y-z
// space of less than the robot's Precision.
func (r *Robot) linear(start, target Pose, a geom.Angle, v geom.Vector, paces chan<- Pace, total float64, quitter <-chan struct{}) {
	defer close(paces)
	wp0 := Pose{
		J: make([]geom.Angle, len(target.J)),
		M: geom.M(start.M...),
		V: geom.V(start.V...),
	}
	copy(wp0.J[:], start.J[:])
	wpN := target

	frac := 0.0
	top := 1.0
	for {
		if geom.Zeroish(frac - top) {
			paces <- Pace{Frac: -1}
			return
		}
		half := 0.5 * (frac + top)
		mid := Pose{
			V: start.V.AddS(target.V.Sub(start.V), half),
		}
		if a != 0 {
			ang := a * geom.Angle(half)
			rot, err := v.RV(ang)
			if err != nil {
				paces <- Pace{Frac: -1}
				return
			}
			mid.M = rot.XM(start.M)
		} else {
			mid.M = start.M
		}
		js := r.Inverse(mid.M, mid.V)

		var err error
		var n int
		if frac < 0.5 {
			n, err = r.Closest(&start, js)
		} else {
			n, err = r.Closest(&target, js)
		}
		if err != nil {
			paces <- Pace{Frac: -1}
			return
		}
		mid.J = js[n]

		var d geom.Angle
		ok := true
		for i := 0; ok && i < len(r.p.J); i++ {
			d += geom.Angle(math.Abs(float64(mid.J[i] - 0.5*(wp0.J[i]+wpN.J[i]))))
			ok = d < r.Precision
		}
		if !ok {
			wpN = mid
			top = half
		} else if geom.Zeroish(top - 1) {
			break
		} else {
			paces <- Pace{
				Frac: total * top,
				J:    wpN.J,
			}
			frac = top
			top = 1
			wp0 = wpN
			wpN = target
		}
	}
	paces <- Pace{
		Frac: total,
		J:    target.J,
	}
}

// align determines the axis of rotation from the starting pose that
// tracks from the start to the target pose through the planned
// motion.
func (r *Robot) align(was *Pose, target Pose) (Pose, geom.Angle, geom.Vector, error) {
	// The motion will perform this rotation as it moves from its
	// initial pose to the target.
	if was == nil {
		was = &r.p
	}

	// duplicate was - it won't be valid for long if it came from
	// the robot since the robot will likely be stepping while we
	// solve more of its path.
	start := Pose{
		J: make([]geom.Angle, len(was.J)),
		M: geom.M(was.M...),
		V: geom.V(was.V...),
	}
	copy(start.J[:], was.J[:])
	var a geom.Angle
	var v geom.Vector
	if !r.p.M.Equals(target.M) {
		// This motion involves a change of basis.
		inv, err := start.M.Inv()
		if err != nil {
			return Pose{}, 0, nil, err
		}
		change := target.M.XM(inv)
		if _, v, a, err = change.Eigen(); err != nil {
			return Pose{}, 0, nil, err
		}
	}
	return start, a, v, nil
}

// Linear returns a channel over which a backgrounded function
// forwards a sequence of paces that cause the robot to move linearly
// between the was pose and the target pose. The function terminates
// and closes the channel to signal early exit if the quitter channel
// is closed. Total is the range over which we report fractional
// progress at each pose.
func (r *Robot) Linear(was *Pose, total float64, target Pose, quitter <-chan struct{}) (<-chan Pace, error) {
	start, a, v, err := r.align(was, target)
	if err != nil {
		return nil, err
	}
	paces := make(chan Pace, 3)
	go r.linear(start, target, a, v, paces, total, quitter)
	return paces, nil
}

// Joined translates a destination pose into a single Pace over the
// returned channel. This is the computationally simplest form of
// trajectory since the whole path is a uniform Joint interpolation.
func (r *Robot) Joined(was *Pose, total float64, target Pose, quitter <-chan struct{}) (<-chan Pace, error) {
	paces := make(chan Pace)
	go func() {
		defer close(paces)
		select {
		case <-quitter:
		case paces <- Pace{Frac: total, J: target.J}:
		}
	}()
	return paces, nil
}

// arc derives a set of paces that can be each linearly interpolated
// to map out an arc between start and target with the path smoothly
// rotating around a point pivot by a total angle of theta. The basis
// of the faceplate of the robot rotates independently around v by a
// total angle a.
func (r *Robot) arc(start, target Pose, a geom.Angle, v geom.Vector, paces chan<- Pace, total float64, axis, pivot geom.Vector, theta geom.Angle, quitter <-chan struct{}) {
	defer close(paces)
	if geom.Zeroish(theta.Rad()) {
		panic("r.arc called with theta zero - should redirect to linear")
	}
	spoke := start.V.Sub(pivot)

	wp0 := Pose{
		J: make([]geom.Angle, len(target.J)),
		M: geom.M(start.M...),
		V: geom.V(start.V...),
	}
	copy(wp0.J[:], start.J[:])
	wpN := target

	frac := 0.0
	top := 1.0
	for {
		if geom.Zeroish(frac - top) {
			paces <- Pace{Frac: -1}
			return
		}
		half := 0.5 * (frac + top)

		rAxis, err := axis.RV(theta * geom.Angle(half))
		if err != nil {
			// Just terminate the paces.
			paces <- Pace{Frac: -1}
			return
		}
		mid := Pose{
			V: pivot.Add(rAxis.XV(spoke)),
		}
		if a != 0 {
			ang := a * geom.Angle(half)
			rot, err := v.RV(ang)
			if err != nil {
				paces <- Pace{Frac: -1}
				return
			}
			mid.M = rot.XM(start.M)
		} else {
			mid.M = start.M
		}
		js := r.Inverse(mid.M, mid.V)

		var n int
		if frac < 0.5 {
			n, err = r.Closest(&start, js)
		} else {
			n, err = r.Closest(&target, js)
		}
		if err != nil {
			// Ran out of solutions.
			paces <- Pace{Frac: -1}
			return
		}
		mid.J = js[n]

		var d geom.Angle
		ok := true
		for i := 0; ok && i < len(r.p.J); i++ {
			d += geom.Angle(math.Abs(float64(mid.J[i] - 0.5*(wp0.J[i]+wpN.J[i]))))
			ok = d < r.Precision
		}
		if !ok {
			wpN = mid
			top = half
		} else if geom.Zeroish(top - 1) {
			break
		} else {
			paces <- Pace{
				Frac: total * top,
				J:    wpN.J,
			}
			frac = top
			top = 1
			wp0 = wpN
			wpN = target
		}
	}
	paces <- Pace{
		Frac: total,
		J:    target.J,
	}
}

// Arc returns a channel over which a set of paces map out a linearly
// interpolatable set of joint poses that map out an arc from the was
// pose to the target pose. The plane of said arc is perpendicular to
// axis, and the angle of the arc is theta.
func (r *Robot) Arc(was *Pose, total float64, target Pose, axis geom.Vector, theta geom.Angle, quitter <-chan struct{}) (<-chan Pace, error) {
	if geom.Zeroish(theta.Rad()) {
		// Treat this as a straight line.
		return r.Linear(was, total, target, quitter)
	}

	start, a, v, err := r.align(was, target)
	if err != nil {
		return nil, err
	}

	d := target.V.Sub(start.V)
	if d.Equals(geom.ZeroV) {
		return nil, geom.ErrNormalNotPossible
	}

	half := theta * 0.5
	radius := d.R() * 0.5 / half.S()
	n, err := axis.Cross(d).Normalize()
	if err != nil {
		return nil, err
	}
	pivot := target.V.Add(start.V).Scale(0.5).AddS(n, radius*half.C())

	paces := make(chan Pace, 3)
	go r.arc(start, target, a, v, paces, total, axis, pivot, theta, quitter)
	return paces, nil
}

// near performs an interpolated path close to three poses. The path,
// parameterized by a progress parameter, can start arbitrarily far
// through the flight >= 0, and terminate at a later point <= 1.
func (r *Robot) near(a, b, c Pose, from, to float64, theta geom.Angle, axis geom.Vector, paces chan<- Pace, quitter <-chan struct{}) {
	defer close(paces)

	rot, err := axis.RV(theta * geom.Angle(from))
	if err != nil {
		paces <- Pace{Frac: -1}
		return
	}

	oMFrom := 1 - from
	wp0 := Pose{
		M: rot.XM(a.M),
		V: a.V.Scale(oMFrom*oMFrom).AddS(b.V, 2*from*oMFrom).AddS(c.V, from*from),
	}
	js := r.Inverse(wp0.M, wp0.V)
	n, err := r.Closest(&b, js)
	if err != nil {
		paces <- Pace{Frac: -1}
		return
	}
	wp0.J = js[n]

	rot, err = axis.RV(theta * geom.Angle(to))
	if err != nil {
		paces <- Pace{Frac: -1}
		return
	}

	oMTo := 1 - to
	wpE := Pose{
		M: rot.XM(a.M),
		V: a.V.Scale(oMTo*oMTo).AddS(b.V, 2*to*oMTo).AddS(c.V, to*to),
	}
	js = r.Inverse(wpE.M, wpE.V)
	n, err = r.Closest(&b, js)
	if err != nil {
		paces <- Pace{Frac: -1}
		return
	}
	wpE.J = js[n]

	wpN := wpE

	frac, top := from, to
	for {
		if geom.Zeroish(frac - top) {
			paces <- Pace{Frac: -1}
			return
		}
		half := 0.5 * (frac + top)

		rot, err := axis.RV(theta * geom.Angle(half))
		if err != nil {
			// Just terminate the paces.
			paces <- Pace{Frac: -1}
			return
		}

		oMP := (1 - half)
		v := a.V.Scale(oMP*oMP).AddS(b.V, 2*half*oMP).AddS(c.V, half*half)
		mid := Pose{
			V: v,
			M: rot.XM(a.M),
		}
		js := r.Inverse(mid.M, mid.V)
		n, err := r.Closest(&b, js)
		if err != nil {
			paces <- Pace{Frac: -1}
			return
		}
		mid.J = js[n]

		var d geom.Angle
		ok := true
		for i := 0; ok && i < len(r.p.J); i++ {
			d += geom.Angle(math.Abs(float64(mid.J[i] - 0.5*(wp0.J[i]+wpN.J[i]))))
			ok = d < r.Precision
		}
		if !ok {
			wpN = mid
			top = half
		} else if geom.Zeroish(top - to) {
			break
		} else {
			paces <- Pace{
				Frac: top,
				J:    wpN.J,
			}
			frac = top
			top = to
			wp0 = wpN
			wpN = wpE
		}
	}
	paces <- Pace{
		Frac: top,
		J:    wpE.J,
	}
}

// glide stitches together three segments to take total time to
// traverse.
func (r *Robot) glide(u, a, x, y, b, d Pose, ang geom.Angle, axis geom.Vector, ps chan<- Pace, fracX, fracY, total float64, quitter <-chan struct{}) {
	defer close(ps)

	p1 := make(chan Pace, 3)
	go r.near(u, a, x, 0.5, 1, 2*ang*geom.Angle(fracX), axis, p1, quitter)
	for p := range p1 {
		if p.Frac == -1 {
			ps <- p
			return
		}
		p.Frac = (p.Frac - 0.5) * 2 * fracX * total
		ps <- p
	}

	p2 := make(chan Pace, 3)
	go r.linear(x, y, ang*geom.Angle(fracY-fracX), axis, p2, 1, quitter)
	for p := range p2 {
		if p.Frac == -1 {
			ps <- p
			return
		}
		p.Frac = (fracX + p.Frac*(fracY-fracX)) * total
		ps <- p
	}

	p3 := make(chan Pace, 3)
	dFrac := 1 - fracY
	go r.near(y, b, d, 0, 0.5, ang*2*geom.Angle(dFrac), axis, p3, quitter)
	for p := range p3 {
		if p.Frac == -1 {
			ps <- p
			return
		}
		p.Frac = (fracY + p.Frac*2*dFrac) * total
		ps <- p
	}
}

// Near determines a multi-segment path based on three consecutive
// target poses (a,b,c) and two nearness thresholds parameterizing how
// much to cut off the corners (rA, rB). The start location, w, need
// not be any of these poses. The result of the generated paces is to
// adopt a pose near to b (at z), with b's M value.
//
// The path goes, w (curve to) x (line to) y (curve to) z. Where xy
// follows the middle segment of ab. If w is on (ab) then the (curve
// to) x is interpreted as linear.
//
// In all cases, the final location, z, is such that stitching
// together calls to Near will yield continuous motion. That is, z is
// computed purely based on a,b,c and not w.
//
// For each rounded corner, we do not use an arc, but a parameterized
// sum of the nearby vectors. The parameterization has w and z at the
// 50% parameterization.
func (r *Robot) Near(w *Pose, total float64, a, b, c Pose, rA, rB float64, quitter <-chan struct{}) (<-chan Pace, error) {
	// Figure out the starting pose and the basis rotation over
	// the full course of the motion between x and z (z.M=b.M). We
	// will preserve this rotation and angular progress through
	// the course of the combined trajectory.
	start, ang, axis, err := r.align(w, b)
	if err != nil {
		return nil, err
	}
	if geom.Zeroish(float64(ang)) {
		axis = geom.Z(1)
	}

	// The theory of cutting a corner from: j->k->l is we follow
	// the trajectory:
	//
	//    k + (1-p)^2(j-k) + p^2(l-k) = j*(1-p)^2 + 2k*p(1-p) + l*p^2
	//
	// at the half way point, p=0.5: p^2=0.25=(1-p)^2, and the point:
	//
	//    h = 0.25*j + 0.5*k + 0.25*l
	//
	// so, given an h, we can deduce the effective ghost of the
	// prior corner as:
	//
	//    j = 4 * (h - 0.5*k - 0.25*l)
	//
	// For the first half corner, we have the half way point of w
	// (our h), a (our k), and we need x (our l). To compute the
	// ghost of a prior corner, we must first determine x. While
	// we are at it, we compute y and the ghost of d (the point on
	// bc we would intercept after point z).

	ab := b.V.Sub(a.V)
	rAB := ab.R()
	hAB := 0.5 * rAB

	fracX := 0.5
	if hAB > rA {
		fracX = rA / rAB
	}
	rot, err := axis.RV(geom.Angle(fracX) * ang)
	if err != nil {
		return nil, err
	}
	x := Pose{
		V: a.V.AddS(ab, fracX),
		M: rot.XM(start.M),
	}
	js := r.Inverse(x.M, x.V)
	n, err := r.Closest(&a, js)
	if err != nil {
		return nil, err
	}
	x.J = js[n]

	fracY := 0.5
	if hAB > rB {
		fracY = 1 - rB/rAB
	}
	rot, err = axis.RV(geom.Angle(fracY) * ang)
	if err != nil {
		return nil, err
	}
	y := Pose{
		V: a.V.AddS(ab, fracY),
		M: rot.XM(start.M),
	}
	js = r.Inverse(y.M, y.V)
	n, err = r.Closest(&b, js)
	if err != nil {
		return nil, err
	}
	y.J = js[n]

	// We have three segments, with two ghost edges:
	//    - u -> (0) w->x (fracX), x->y (fracY), y->z (1) -> d
	//
	// u is our j for the first corner, recall we are at w (aka
	// start) already, so u is backwards and only useful as a
	// reference point and a reference orientation. d is the
	// destination we will only get half way to.
	rot, err = axis.RV(-geom.Angle(fracX) * ang)
	if err != nil {
		return nil, err
	}
	u := Pose{
		V: start.V.AddS(a.V, -0.5).AddS(x.V, -0.25).Scale(4),
		M: rot.XM(start.M),
	}
	js = r.Inverse(u.M, u.V)
	n, err = r.Closest(&a, js)
	if err != nil {
		return nil, err
	}
	u.J = js[n]

	dFrac := 1 - fracY
	rot, err = axis.RV(geom.Angle(1+dFrac) * ang)
	if err != nil {
		return nil, err
	}

	bc := c.V.Sub(b.V)
	rBC := bc.R()
	hBC := 0.5 * rBC
	fracD := 0.5
	if hBC > rB {
		fracD = rB / rBC
	}
	d := Pose{
		V: b.V.AddS(bc, fracD),
		M: rot.XM(start.M),
	}
	js = r.Inverse(d.M, d.V)
	n, err = r.Closest(&b, js)
	if err != nil {
		return nil, err
	}
	d.J = js[n]

	ps := make(chan Pace, 3)
	go r.glide(u, a, x, y, b, d, ang, axis, ps, fracX, fracY, total, quitter)
	return ps, nil
}
