package org.sciencework.core;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.ode.FirstOrderIntegrator;

/**
 * CVRA means Complex Value and Real Argument
 * It is a complex number adapter for {@link FirstOrderIntegrator}
 * @author akutlin
 *
 */

public interface CVRAFirstOrderIntegrator {
	public double integrate(CVRAFirstOrderDifferentialEquations equations, double t0, Complex[] y0, double t, Complex[] y);
	public void setFirstOrderIntegrator( FirstOrderIntegrator integrator );
}
