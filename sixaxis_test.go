package sixaxis

import (
	"testing"

	"zappem.net/pub/math/geom"
)

func TestForward(t *testing.T) {
	js := []Joint{
		{Min: geom.Degrees(-170), Max: geom.Degrees(170), Rot: geom.RZ, Offset: 10},
		{Min: geom.Degrees(-120), Max: geom.Degrees(120), Rot: geom.RX, Offset: 50},
		{Min: geom.Degrees(-120), Max: geom.Degrees(120), Rot: geom.RX, Offset: 30},
		{Min: geom.Degrees(-170), Max: geom.Degrees(170), Rot: geom.RZ, Offset: 20},
		{Min: geom.Degrees(-120), Max: geom.Degrees(120), Rot: geom.RX, Offset: 10},
		{Min: geom.Degrees(-360), Max: geom.Degrees(360), Rot: geom.RZ, Offset: 0},
	}
	reach := 0.0
	for _, j := range js {
		reach += j.Offset
	}

	r, err := NewRobot(js...)
	if err != nil {
		t.Fatalf("failed to define new robot: %v", err)
	}
	m, v := r.MV()
	if !m.Equals(geom.I) {
		t.Errorf("non-identity for un-posed robot: %v", m)
	}
	if want := geom.Z(reach); !v.Equals(want) {
		t.Errorf("fwd kinematics location for un-posed robot: got=%v want=%v", v, want)
	}

	// Pose the robot in some non-trivial poses and validate that the best
	// ik solution matches the arranged pose.
	for i := range js {
		if err := r.SetJ(i, geom.Degrees(11+2*float64(i))); err != nil {
			t.Errorf("failed to set joint[%d] to %d deg", i, 13)
		}
		m, v := r.MV()
		jjs := r.Inverse(m, v)
		if len(jjs) == 0 {
			t.Fatalf("unable to solve J[%d] inverse", i)
		}
		if best, err := r.Closest(nil, jjs); err != nil {
			t.Errorf("inverse yielded %d solutions: %v", len(jjs), err)
		} else {
			for k, a := range jjs[best] {
				if !geom.Zeroish(float64(a - r.p.J[k])) {
					t.Errorf("[%d] ik solution: %v != %v", i, jjs[best], r.p.J)
					break
				}
			}
		}
	}
}
