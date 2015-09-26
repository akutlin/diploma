package org.sciencework.core;

import org.apache.commons.math3.exception.DimensionMismatchException;

/**
 * CVRA means Complex Value and Real Argument
 * It is a complex number adapter for {@link UnivariateDifferentiableFunction}
 * @author akutlin
 *
 */
public interface CVRAUnivariateDifferentiableFunction extends CVRAUnivariableFunction {
	ComplexDerivativeStructure value(ComplexDerivativeStructure t)
            throws DimensionMismatchException;
}
