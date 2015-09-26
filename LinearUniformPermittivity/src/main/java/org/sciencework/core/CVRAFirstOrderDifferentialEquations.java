package org.sciencework.core;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

/**
 * CVRA means Complex Value and Real Argument
 * It is a complex number adapter for {@link FirstOrderDifferentialEquations}
 * @author akutlin
 *
 */

public interface CVRAFirstOrderDifferentialEquations {
	
	public FirstOrderDifferentialEquations getRealEquations();
	public FirstOrderDifferentialEquations getImaginaryEquations();
	
	/**
	 * Get the dimension of the problem.
	 * @return
	 */
	public int getDimension();
}
