package org.sciencework;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.junit.Test;

public class TETest {
	
	DormandPrince853Integrator integrator = new  DormandPrince853Integrator(1E-17, 1E-3, 1E-10, 1E-10);
	
	@Test
	public void vacuum() {
		
		Permittivity vacuum = new Permittivity() {

			public double getReal(double x) {
				return 1;
			}

			public double getImg(double x) {
				return 0;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(1);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(0);
			}
			
		};
		
		double teta = 0;
		double w = 1E11;
		
		TE wave = new TE(teta, vacuum, w);
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = -100;
		double t1 = 100;
		integrator.integrate(wave, t0, y, t1, y);
		
		RealVector S0 = wave.getPowerFlow(Ey0, Bz0);
		double S0N = S0.getNorm();
		
		RealVector S1 = wave.getPowerFlow( new Complex( y[0], y[1] ), new Complex( y[2], y[3] ) );
		double S1N = S1.getNorm();
		
		double delta = ( S1N - S0N )/S0N;
		
		assertEquals( 0, delta, 0.001);
	}
	
	@Test
	public void vacuumWithAngle() {
		
		Permittivity vacuum = new Permittivity() {

			public double getReal(double x) {
				return 1;
			}

			public double getImg(double x) {
				return 0;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(1);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(0);
			}
			
		};
		
		double teta = 1;
		double w = 1E11;
		
		TE wave = new TE(teta, vacuum, w);
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = -100;
		double t1 = 100;
		integrator.integrate(wave, t0, y, t1, y);
		
		RealVector S0 = wave.getPowerFlow(Ey0, Bz0);
		double S0N = S0.getNorm();
		
		RealVector S1 = wave.getPowerFlow( new Complex( y[0], y[1] ), new Complex( y[2], y[3] ) );
		double S1N = S1.getNorm();
		
		double delta = ( S1N - S0N )/S0N;
		
		assertEquals( 0, delta, 0.001);
	}
	
	@Test
	public void stepFunctionPermittivity() {
		
		Permittivity vacuum = new Permittivity() {

			public double getReal(double x) {
				return x > 0 ? 2 : 1;
			}

			public double getImg(double x) {
				return 0;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(1);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(0);
			}
			
		};
		
		double teta = 0;
		double w = 1E11;
		
		TE wave = new TE(teta, vacuum, w);
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = -100;
		double t1 = 100;
		integrator.integrate(wave, t0, y, t1, y);
		
		RealVector S0 = wave.getPowerFlow(Ey0, Bz0);
		double S0N = S0.getNorm();
		
		RealVector S1 = wave.getPowerFlow( new Complex( y[0], y[1] ), new Complex( y[2], y[3] ) );
		double S1N = S1.getNorm();
		
		double delta = ( S0N - S1N )/S1N;
		
		assertEquals( 0, delta, 0.001);
		assertEquals( 1, S1N, 0.001);

	}
	
	@Test
	public void stepFunctionPermittivityAngle() {
		
		Permittivity vacuum = new Permittivity() {

			public double getReal(double x) {
				return x > 0 ? 2 : 1;
			}

			public double getImg(double x) {
				return 0;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(1);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(0);
			}
			
		};
		
		double teta = 1;
		double w = 1E11;
		
		TE wave = new TE(teta, vacuum, w);
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = 100;
		double t1 = -100;
		integrator.integrate(wave, t0, y, t1, y);
		
		RealVector S0 = wave.getPowerFlow(Ey0, Bz0);
		double S0N = S0.getNorm();
		
		RealVector S1 = wave.getPowerFlow( new Complex( y[0], y[1] ), new Complex( y[2], y[3] ) );
		double S1N = S1.getNorm();
		
		double delta = ( S0N - S1N )/S1N;
		
		assertEquals( 0, delta, 0.001);
		assertEquals( 1, S1N, 0.001);
	}
	
	@Test
	public void dissipativeVacuum() {
		
		final double real = 1;
		final double imag = 0.0001;
		double k = 10;
		
		Permittivity vacuum = new Permittivity() {

			public double getReal(double x) {
				return real;
			}

			public double getImg(double x) {
				return imag;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(real);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(imag);
			}
			
		};
		
		double teta = 0;
		
		TE wave = new TE(teta, vacuum, k * TE.c );
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = -100;
		double t1 = 100;
		integrator.integrate(wave, t0, y, t1, y);
		
		RealVector S0 = wave.getPowerFlow(Ey0, Bz0);
		double S0N = S0.getNorm();
		
		RealVector S1 = wave.getPowerFlow( new Complex( y[0], y[1] ), new Complex( y[2], y[3] ) );
		double S1N = S1.getNorm();
		
		double delta = ( S1N - S0N )/S0N;
		
		Complex amp = new Complex( real - Math.pow( Math.sin(teta), 2), imag ).pow(0.5).multiply( new Complex( 0, k*( t1 - t0 )) ).exp();
		double theoretical = Math.pow( amp.abs(), 2) - 1;
		
		assertEquals( theoretical, delta, 0.001);
		assertTrue( delta < 0);
	}
	
	@Test
	public void antidissipativeVacuum() {
		
		final double real = 1;
		final double imag = -0.0001;
		double k = 10;
		
		Permittivity vacuum = new Permittivity() {

			public double getReal(double x) {
				return real;
			}

			public double getImg(double x) {
				return imag;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(real);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(imag);
			}
			
		};
		
		double teta = 0;
		
		TE wave = new TE(teta, vacuum, k * TE.c );
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = -100;
		double t1 = 100;
		integrator.integrate(wave, t0, y, t1, y);
		
		RealVector S0 = wave.getPowerFlow(Ey0, Bz0);
		double S0N = S0.getNorm();
		
		RealVector S1 = wave.getPowerFlow( new Complex( y[0], y[1] ), new Complex( y[2], y[3] ) );
		double S1N = S1.getNorm();
		
		double delta = ( S1N - S0N )/S0N;
		
		Complex amp = new Complex( real - Math.pow( Math.sin(teta), 2), imag ).pow(0.5).multiply( new Complex( 0, k*( t1 - t0 )) ).exp();
		double theoretical = Math.pow( amp.abs(), 2) - 1;
		
		assertEquals( theoretical, delta, 0.001);
		assertTrue( delta > 0);
	}
	
	@Test
	public void dissipativeVacuumAngle() {
		
		final double real = 1;
		final double imag = 0.0001;
		double k = 10;
		
		Permittivity vacuum = new Permittivity() {

			public double getReal(double x) {
				return real;
			}

			public double getImg(double x) {
				return imag;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(real);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(imag);
			}
			
		};
		
		double teta = 1;
		
		TE wave = new TE(teta, vacuum, k * TE.c );
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = -100;
		double t1 = 100;
		integrator.integrate(wave, t0, y, t1, y);
		
		RealVector S0 = wave.getPowerFlow(Ey0, Bz0);
		double S0N = S0.getNorm();
		
		RealVector S1 = wave.getPowerFlow( new Complex( y[0], y[1] ), new Complex( y[2], y[3] ) );
		double S1N = S1.getNorm();
		
		double delta = ( S1N - S0N )/S0N;
		
		Complex amp = new Complex( real - Math.pow( Math.sin(teta), 2), imag ).pow(0.5).multiply( new Complex( 0, k*( t1 - t0 )) ).exp();
		double theoretical = Math.pow( amp.abs(), 2) - 1;
		
		assertEquals( theoretical, delta, 0.001);
	}
	
	@Test
	public void solid() {
		
		Permittivity solid = new Permittivity() {

			public double getReal(double x) {
				return Math.abs(x) < 10 * Math.PI ? (1 + Math.cos(x))/2 : 1;
			}

			public double getImg(double x) {
				return 0;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(1);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(0);
			}
			
		};
		
		double teta = 0;
		double w = 1E11;
		
		TE wave = new TE(teta, solid, w);
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = 100;
		double t1 = -100;
		integrator.integrate(wave, t0, y, t1, y);
		
		RealVector S0 = wave.getPowerFlow(Ey0, Bz0);
		double S0N = S0.getNorm();
		
		RealVector S1 = wave.getPowerFlow( new Complex( y[0], y[1] ), new Complex( y[2], y[3] ) );
		double S1N = S1.getNorm();
		
		double delta = ( S0N - S1N )/S1N;
		
		assertEquals( 0, delta, 0.001);
	}
	
	@Test
	public void solidInverce() {
		
		Permittivity solid = new Permittivity() {

			public double getReal(double x) {
				return Math.abs(x) < 10 * Math.PI ? (1 + Math.cos(x))/2 : 1;
			}

			public double getImg(double x) {
				return 0;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(1);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(0);
			}
			
		};
		
		double teta = 0;
		double w = 1E11;
		
		TE wave = new TE(teta, solid, w);
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		RealVector initial = MatrixUtils.createRealVector(y);
		
		double t0 = 100;
		double t1 = -100;
		integrator.integrate(wave, t0, y, t1, y);
		integrator.integrate(wave, t1, y, t0, y);
		
		RealVector finaly = MatrixUtils.createRealVector(y);
		RealVector delta = finaly.subtract(initial);
		
		assertTrue( delta.getNorm() + " > 0.001", delta.getNorm() < 0.001 );
	}
	
	@Test
	public void solidNegative() {
		
		final double size = 10;
		
		Permittivity solid = new Permittivity() {

			public double getReal(double x) {
				return Math.abs(x) < size ? -1 : 1;
			}

			public double getImg(double x) {
				return 0;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(1);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(0);
			}
			
		};
		
		double teta = 0;
		double k = 1;
		
		TE wave = new TE(teta, solid, k * TE.c);
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = 20;
		double t1 = -20;
		integrator.integrate(wave, t0, y, t1, y);
		
		RealVector S0 = wave.getPowerFlow(Ey0, Bz0);
		double S0N = S0.getNorm();
		
		RealVector S1 = wave.getPowerFlow( new Complex( y[0], y[1] ), new Complex( y[2], y[3] ) );
		double S1N = S1.getNorm();
		
		double delta = ( S0N - S1N )/S1N;
		
		assertEquals( 0, delta, 0.001);
	}
	
	@Test
	public void solidNegativeInverse() {
		
		final double size = 10;
		
		Permittivity solid = new Permittivity() {

			public double getReal(double x) {
				return Math.abs(x) < size ? -1 : 1;
			}

			public double getImg(double x) {
				return 0;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(1);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(0);
			}
			
		};
		
		double teta = 0;
		double k = 1;
		
		TE wave = new TE(teta, solid, k * TE.c);
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		RealVector initial = MatrixUtils.createRealVector(y);
		
		double t0 = 20;
		double t1 = -20;
		integrator.integrate(wave, t0, y, t1, y);
		integrator.integrate(wave, t1, y, t0, y);
		
		RealVector finaly = MatrixUtils.createRealVector(y);
		RealVector delta = finaly.subtract(initial);
		
		assertTrue( delta.getNorm() + " > 0.001", delta.getNorm() < 0.001 );
		
	}
	
	@Test
	public void solidNegativeWithoutField() {
		
		final double size = 10;
		
		Permittivity solid = new Permittivity() {

			public double getReal(double x) {
				return Math.abs(x) < size ? -1 : 1;
			}

			public double getImg(double x) {
				return 0;
			}

			public DerivativeStructure getRealDerivative(DerivativeStructure x) {
				return x.createConstant(1);
			}

			public DerivativeStructure getImgDerivative(DerivativeStructure x) {
				return x.createConstant(0);
			}
			
		};
		
		double teta = 0;
		double k = 1;
		
		TE wave = new TE(teta, solid, k * TE.c);
		
		Complex Ey0 = new Complex(0,0);
		Complex Bz0 = new Complex(0 * Math.cos(teta),0);
						
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		RealVector initial = MatrixUtils.createRealVector(y);
		
		double t0 = 20;
		double t1 = -20;
		integrator.integrate(wave, t0, y, t1, y);
		
		RealVector finaly = MatrixUtils.createRealVector(y);
		RealVector delta = finaly.subtract(initial);
		
		assertEquals( 0, finaly.getNorm(), 0);
		assertEquals( 0, delta.getNorm(), 0);

	}
}
