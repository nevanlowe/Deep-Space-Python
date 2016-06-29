import math

G = 6.67408e-11
g0 = 9.81
KERBINMASS = 5.2915793e+22


def deltaV(isp, fuel):
    """Return the DeltaV of a rocket's stage.

    isp -- Specific Impulse of the rocket engine (in seconds)
    fuel -- Mass fraction of that stage
    """
    # Tsiolkovsky Rocket Equation
    deltav = isp * g0 * math.log(fuel)
    return deltav


def wetDryRatio(isp, dv):
    """Return the mass fraction needed to reach a given DeltaV.

    Note that as of KSP version 1.1.2, a ratio > 9 is not possible using stock
    fuel tanks. If you get an answer higher than 9, then you are SOL.

    isp -- Specific Impulse of the rocket engine (in seconds)
    dv -- DeltaV you want to reach
    """
    # Tsiolkovsky Rocket Equation rearranged to solve for mass fraction.
    ratio = math.exp(dv / (isp * g0))
    return ratio


def semiMajorAxis(peri, apo, r):
    """Return the semi-major axis of an orbit.

    peri -- Periapsis of your orbit
    apo -- Apoapsis of your orbit
    r -- Radius of central body
    """
    # Adds the central body's radius to the two altitudes because KSP's
    # altitude meter does not include the central body's radius.
    sma = (peri + apo) / 2 + r
    return sma


def semiMajorApo(sma, peri, r):
    """Return the apoapsis needed to reach a given semi-major axis.

    sma -- Semi-major axis
    peri -- Periapsis of your orbit
    r -- Radius of central body
    """
    # semiMajorAxis rearranged to solve for apoapsis
    apo = 2 * (sma - r) - peri
    return apo


def semiMajorPeri(sma, apo, r):
    """Return the apoapsis needed to reach a given semi-major axis.

    sma -- Semi-major axis
    apo -- Apoapsis of your orbit
    r -- Radius of central body
    """
    # semiMajorAxis rearranged to solve for periapsis
    peri = 2 * (sma - r) - apo
    return peri


def orbitalPeriod(sma, M):
    """Return the sidereal period of an orbit.

    sma -- Semi-major axis of orbiting body
    M -- Mass of central body
    """
    # Kepler's Third Law arranged to solve for orbital period
    T = 2 * math.pi * (sma ** 3 / (G * M)) ** 0.5
    return T


def semiMajorPeriod(T, M):
    """Return the semi-major axis needed for a given orbital period

    T -- Sidereal orbital period
    M -- Mass of central body
    """
    # Kepler's Third Law arranged to solve for semi-major axis
    sma = (G * M * T ** 2 / (4 * math.pi ** 2)) ** (1 / 3)
    return sma


def eccentricity(peri, apo, r):
    """Return the eccentricity of an elliptical orbit.

    peri -- Periapsis of your orbit
    apo -- Apoapsis of your orbit
    r -- Radius of central body
    """
    # Adds the central body's radius to the two altitudes because KSP's
    # altitude meter does not include the central body's radius.
    maximum = apo + r
    minimum = peri + r
    # Elliptical orbit eccentricity formula
    e = (maximum - minimum) / (maximum + minimum)
    return e


def hyperEccentricity(peri, r, v, M):
    """Return the eccentricity of a hyperbolic orbit.

    The eccentricity() function can't be used for hyperbolic orbits because
    hyperbolic orbits lack an apoapsis. This function uses other orbital
    parameters instead, removing the need for an apoapsis. This can also be
    used for parabolic orbits, but I don't know why you would ever need to
    calculate the eccentricity of a parabolic orbit ever. Seriously. If you
    ever need to calculate the eccentricity of a parabolic orbit, you are doing
    something horribly horribly wrong and you should stop to reconsider your
    life and/or conic sections.

    peri -- Periapsis of hyperbolic orbit
    r -- Radius of central body
    v -- Velocity at periapsis
    M -- Mass of central body
    """
    # Adds the central body's radius to the periapsis because KSP's
    # altitude meter does not include the central body's radius.
    h = peri + r
    # Vis-viva rearranged to solve for semi-major axis
    sma = h * G * M / (2 * G * M - h * v ** 2)
    # Eccentricity of a hyperbola formula
    e = (sma - h) / sma
    return e


def hohmannAngle(h1, h2, r):
    """Return the phase angle between you and the target in a Hohmann transfer.

    This is necessary for if you want to actually line up with the target at
    the end and have an encounter.

    h1 -- Initial orbit height
    h2 -- Final orbit height
    r -- Radius of central body
    """
    # Via Kepler's Third Law, this is number of orbits the target body
    # completes during half of the elliptical transfer orbit.
    pt = 0.5 * ((h1 + h2 + 2 * r) / (2 * r + 2 * h2)) ** 1.5
    # Ignoring all digits to the left of the decimal (thus disregarding
    # complete orbits), this is the number of degrees that the target body will
    # travel during half of the elliptical transfer orbit.
    ft = pt % 1
    sweep = ft * 360
    # The point in your orbit that you should burn so that you will line up
    # with the target at the apoapsis of the elliptical transfer orbit.
    phase = 180 - sweep
    return phase


def hohmannVelocity(h1, h2, r, M, singlevalue=False):
    """Return necessary DeltaV required for a Hohmann Transfer.

    Can either return the value for both burns as a tuple (more useful for
    orbital maneuvering) or as the sum of the burns as a single value (more
    convenient for DeltaV budget analysis).

    h1 -- Initial orbit height
    h2 -- Final orbit height
    r -- Radius of central body
    M -- Mass of central body
    singlevalue -- True returns the sum of the DeltaV needed, False (default)
    returns a tuple
    """
    # Adds the central body's radius to the two altitudes because KSP's
    # altitude meter does not include the central body's radius.
    r1 = h1 + r
    r2 = h2 + r
    # This calculates the DeltaV for the first burn (the one that creates an
    # elliptical orbit that starts the transfer), and the second burn (the one
    # that circularizes your orbit at the apoaspis of the elliptical orbit).
    dv1 = (G * M / r1) ** 0.5 * ((2 * r2 / (r1 + r2)) ** 0.5 - 1)
    dv2 = (G * M / r2) ** 0.5 * (1 - (2 * r1 / (r1 + r2)) ** 0.5)
    if singlevalue:
        return dv1 + dv2
    else:
        return dv1, dv2


def escapeSurface(M, r):
    """Return escape velocity from the surface of a body.

    M -- Mass of central body
    r -- Radius of central body
    """
    # Degenerate Vis-Viva with infinite SMA (or equivalently, 1/a = 0).
    ve = (2 * G * M / r) ** 0.5
    return ve


def orbitalVelocity(M, r, h, a):
    """Return orbital speed via the Vis-Viva equation.

    M -- Mass of central body
    r -- Radius of central body
    h -- Current height above surface
    a -- Semi-major axis of orbit
    """
    # Adds the central body's radius to the current height because KSP's
    # altitude meter does not include the central body's radius.
    R = h + r
    # Vis-Viva for a given orbital height
    v = (G * M * (2 / R - 1 / a)) ** 0.5
    return v


def escapeOrbit(M, r, h, a):
    """Return the DeltaV required to escape from orbit.

    M -- Mass of central body
    r -- Radius of central body
    h -- Current height above surface
    a -- Semi-major axis of orbit
    """
    # Adds the central body's radius to the current height because KSP's
    # altitude meter does not include the central body's radius.
    R = h + r
    # Escape velocity from a given height minus Vis-Viva from the same height.
    ve = (2 * G * M / R) ** 0.5 - orbitalVelocity(M, r, h, a)
    return ve


def ejectionVelocity(hev, M, peri, r, delta=False, sma=None):
    """Return the velocity at periapsis for a given Hyperbolic Excess Velocity.

    It can also alternatively return the DeltaV required to achieve that
    ejection velocity. If you want the necessary DeltaV, then the parking orbit
    assumed to be circular unless you provide a semi-major axis.

    Note that this solution is only approximate, as your craft will leave the
    host's sphere of influence before it reaches hyperbolic excess velocity.
    This gives a slightly larger velocity than actually necessary, but this is
    not a huge deal as having an excess deltaV "safety net" is usually
    necessary anyway because of things like atmospheres and sub-optimal
    launches. The timing of the encounter is a problem and will change the
    phase angle, but spheres of influence are fairly large and you can always
    adjust your final orbit.

    hev -- Desired hyperbolic excess velocity
    M -- Mass of central body
    peri -- Periapsis of hyperbolic orbit
    r -- Radius of central body
    delta -- False returns orbital velocity at periapsis, true returns the
    necessary DeltaV
    sma -- Semi-major axis of parking orbit
    """
    # Adds the central body's radius to the periapsis because KSP's
    # altitude meter does not include the central body's radius.
    h = peri + r
    # Hyperbolic excess velocity rearranged to solve for hyperbolic semi-major
    # axis.
    hsma = -G * M / hev ** 2
    # Vis-viva for our hyperbolic orbit at periapsis
    v = ((2 / h - 1 / hsma) * G * M) ** 0.5
    if delta:
        if sma == None:
            sma = h
        # Previous equation minus vis-viva for the parking orbit.
        delta = v - orbitalVelocity(M, r, peri, sma)
        return delta
    else:
        return v


def ejectionVelocity2(ev, SoI, M, peri, r, delta=False, sma=None):
    """Return the velocity at periapsis for a given Excess Velocity.

    This is an attempt to make a more accurate form of the ejectionVelocity
    function by using desired velocity at SoI change, not at infinity like the
    previous function

    ev -- Desired excess velocity leaving the planet
    SoI -- Sphere of Influence
    M -- Mass of central body
    peri -- Periapsis of hyperbolic orbit
    r -- Radius of central body
    delta -- False returns orbital velocity at periapsis, true returns the
    necessary DeltaV
    sma -- Semi-major axis of parking orbit
    """
    # Adds the central body's radius to the periapsis because KSP's
    # altitude meter does not include the central body's radius.
    h = peri + r
    hsma = 1 / (2 / SoI - ev ** 2 / (G * M))
    # Vis-viva for our hyperbolic orbit at periapsis
    v = ((2 / h - 1 / hsma) * G * M) ** 0.5
    if delta:
        if sma is None:
            sma = h
        # Previous equation minus vis-viva for the parking orbit.
        delta = v - orbitalVelocity(M, r, peri, sma)
        return delta
    else:
        return v


def ejectionAngle(hev, M, peri, r):
    """Return the angle required for your HEV to align with the planet's orbit.

    Doing an ejection burn is more efficient than doing a standard Hohmann
    maneuver (because of the Oberth Effect or something), but it's more
    complex. You still need a phase angle, but you also need an ejection angle
    in order for your hyperbolic excess velocity to be parallel with your
    planet's prograde or retrograde vectors. The angle is between the
    transverse axis (basically a line between you and the center of your
    planet) and the planet's prograde or retrograde vector.

    Note that this solution is only approximate, as your craft will leave the
    host's sphere of influence before it reaches hyperbolic excess velocity.
    So the angle is a slight overestimation, and an optimal burn would have a
    slightly smaller angle. This is a larger problem than the velocity because
    an incorrect angle has no advantages, while excess velocity has benefits.
    """
    # Adds the central body's radius to the periapsis because KSP's
    # altitude meter does not include the central body's radius.
    h = peri + r
    # Characteristic Energy (C3) rearranged to solve for hyperbolic semi-major
    # axis. Hyperbolic excess velocity squared is equal to C3
    hsma = -G * M / hev ** 2
    # The hyperbolic parameters a is the distance between the center of the
    # hyperbola and a vertex, b is the vertical distance between the vertex and
    # asymptote, and c is the distance between the center of the hyperbola and
    # a focal point. The three form a right triangle, and c is just a + the
    # distance to the focal point, which is just the periapsis. So the inverse
    # cosine of a/c can be used to find the angle between the transverse axis
    # and the asymptotes (in radians).
    asymptoticAngle = math.acos(hsma / (hsma - h))
    # The angle of ejection and asymptotic angle are supplementary, so we can
    # subtract one from 180 to get the other. We still need to convert to
    # degrees.
    ejectAngle = 180 - math.degrees(asymptoticAngle)
    return ejectAngle


def ejectionAngle2(ev, SoI, M, peri, r):
    """Return the angle required for your EV to align with the planet's orbit.

    This is an incomplete attempt at fixing the inaccuracies in ejectionAngle()
    """
    h = peri + r
    a = 1 / (2 / SoI - ev ** 2 / (G * M))
    c = abs(a - h)
    e = c/-a
    trueAnom = math.acos((a*(1-e**2)-SoI)/(SoI*e))
    fp = ((2*c)**2 + SoI**2 - 2*2*c*SoI*math.cos(trueAnom))**0.5
    A = math.asin(math.sin(trueAnom)*2*c/fp)
    angle = math.pi-trueAnom-A/2
    ejectAngle=180-math.degrees(angle)
    return fp, SoI-2*a, math.degrees(A), ejectAngle


if __name__ == '__main__':
    pass
else:
    print('Hullo, Scott Manley here')