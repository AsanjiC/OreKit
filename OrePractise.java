package org.orekit;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Locale;

import org.hipparchus.analysis.ParametricUnivariateFunction;
import org.hipparchus.fitting.AbstractCurveFitter;
import org.hipparchus.fitting.PolynomialCurveFitter;
import org.hipparchus.fitting.WeightedObservedPoint;
import org.hipparchus.linear.DiagonalMatrix;
import org.hipparchus.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresBuilder;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem;
import org.hipparchus.util.FastMath;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.forces.ForceModel;
import org.orekit.forces.drag.DragForce;
import org.orekit.forces.drag.IsotropicDrag;
import org.orekit.forces.drag.atmosphere.Atmosphere;
import org.orekit.forces.drag.atmosphere.SimpleExponentialAtmosphere;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;

public class OrePractise {
	 
    public static void main(String[] args) throws IOException {
        try {

        		//Set HDD aS main Directory
            File home       = new File(System.getProperty("user.home"));
            File file = new File("OrekitData(Time).txt");
            if(!file.exists())
            {
            	file.createNewFile();
            }
            //Search for OREkit folder in the main directory set above
            File orekitData = new File(home, "orekit-data");
            if (!orekitData.exists()) {
                System.err.format(Locale.US, "Failed to find %s folder%n",
                                  orekitData.getAbsolutePath());
                System.err.format(Locale.US, "You need to download %s from the %s page and unzip it in %s for this tutorial to work%n",
                                  "orekit-data.zip", "https://www.orekit.org/forge/projects/orekit/files",
                                  home.getAbsolutePath());
                System.exit(1);
            }
            
            PrintStream out = new PrintStream(file);
            DataProvidersManager manager = DataProvidersManager.getInstance();
            manager.addProvider(new DirectoryCrawler(orekitData));

            
            // Inertial frame 
            Frame inertialFrame = FramesFactory.getEME2000();

            // Initial date in UTC time scale
            TimeScale utc = TimeScalesFactory.getUTC();
            AbsoluteDate initialDate = new AbsoluteDate(2017, 12, 20, 23, 30, 00.000, utc);

            // gravitation coefficient
            double mu =  3.986004415e+14;
            
            // Initial orbit parameters
            double a = 24396159; // semi major axis in meters
            double e = 0.72831215; // eccentricity
            double i = FastMath.toRadians(7); // inclination
            double omega = FastMath.toRadians(180); // perigee argument
            double raan = FastMath.toRadians(261); // right ascension of ascending node
            double lM = 0; // mean anomaly
            	
            // Orbit construction as Keplerian
            Orbit initialOrbit = new KeplerianOrbit(a, e, i, omega, raan, lM, PositionAngle.MEAN,
                                                    inertialFrame, initialDate, mu);
            //initial State Definition
            SpacecraftState initialState = new SpacecraftState(initialOrbit);
            
            // Step integrator 
            final double minStep = 0.001;
            final double maxstep = 2000.0;
            final double positionTolerance = 10.0;
            final OrbitType propagationType = OrbitType.KEPLERIAN;
            final double[][] tolerances =
                    NumericalPropagator.tolerances(positionTolerance, initialOrbit, propagationType);
            AdaptiveStepsizeIntegrator integrator =
                    new DormandPrince853Integrator(minStep, maxstep, tolerances[0], tolerances[1]);

            // Propagator
            NumericalPropagator propagator = new NumericalPropagator(integrator);
            propagator.setOrbitType(propagationType);

            // Force Model
            final double ae = Constants.GRS80_EARTH_EQUATORIAL_RADIUS;
            OneAxisEllipsoid earth = new OneAxisEllipsoid(ae, Constants.WGS84_EARTH_FLATTENING, FramesFactory.getEME2000());
            earth.setAngularThreshold(1.e-6);
            Atmosphere atm = new SimpleExponentialAtmosphere(earth, 4.e-13, 300000.0, 60000.0);
            final double cd = 2.0;
            final double sf = 5.0;
            ForceModel dragNUM = new DragForce(atm, new IsotropicDrag(sf, cd));

            // Add force model to the propagator
            propagator.addForceModel(dragNUM);

            // Set up initial state in the propagator
            propagator.setInitialState(initialState);

            // Set up operating mode for the propagator as master mode
            // with fixed step(seconds) and specialized step handler
            propagator.setMasterMode(1., new stepHandler());

            // Extrapolate from the initial to the final date
            SpacecraftState finalState = propagator.propagate(initialDate.shiftedBy(7.884e+6));
            /*Updated to Propagate for 1 Month
            SpacecraftState finalState = propagator.propagate(new AbsoluteDate(initialDate, 2629000.));
            */
            KeplerianOrbit o = (KeplerianOrbit) OrbitType.KEPLERIAN.convertType(finalState.getOrbit());
            System.out.print(o.initPVCoordinates().getPosition()+"\n");
            
            //These NextTwo Lines is How I Convert the coordinates to get Long. Lat & Altitude 
            //FieldGeodeticPoint<DerivativeStructure> gp = earth.transform(o.initPVCoordinates(),inertialFrame,initialDate.shiftedBy(10.));
            //System.out.print(gp.toString() +"\n");
            
            // estimator
            final BatchLSEstimator estimator = createEstimator(parser, propagator);
            
            
            System.out.format(Locale.US, "Final state:%n%s %12.3f %10.8f %10.6f %10.6f %10.6f %10.6f%n",
                              finalState.getDate(),
                              o.getA(), o.getE(),
                              FastMath.toDegrees(o.getI()),
                              FastMath.toDegrees(o.getPerigeeArgument()),
                              FastMath.toDegrees(o.getRightAscensionOfAscendingNode()),
                              FastMath.toDegrees(o.getTrueAnomaly()));
     	
         out.close();
        } catch (OrekitException oe) {
            System.err.println(oe.getMessage());
        }
    }
    
    
    // Step Handler Class used to print on the output stream at the given step.
    private static class stepHandler implements OrekitFixedStepHandler {
 
    	
        private stepHandler() {
            //private constructor
        }
        
        public  void init(final SpacecraftState s0, final AbsoluteDate t, final double step) {
        		System.out.println("          p	           	\tv" +
                               "                     a");
        }
        public  void handleStep(SpacecraftState currentState, boolean isLast){
        	
        	KeplerianOrbit o = (KeplerianOrbit) OrbitType.KEPLERIAN.convertType(currentState.getOrbit());
        	System.out.print(o.getDate() + " ");
            System.out.print(o.initPVCoordinates().getPosition()+"\n"); 
            
            
            /*
            System.out.format(Locale.US, "%s %12.3f %10.8f %10.6f %10.6f %10.6f %10.6f%n",
                              currentState.getDate(),
                              o.getA(), o.getE(),
                              FastMath.toDegrees(o.getI()),
                              FastMath.toDegrees(o.getPerigeeArgument()),
                              FastMath.toDegrees(o.getRightAscensionOfAscendingNode()),
                              FastMath.toDegrees(o.getTrueAnomaly()));
            */
            if (isLast) {
                System.out.println("this was the last step ");
                System.out.println();
            }

        }
        	
    }
    
    public class SecularAndHarmonic {

        /** Degree of polynomial secular part. */
        private final int secularDegree;

        /** Pulsations of harmonic part. */
        private final double[] pulsations;

        /** Reference date for the model. */
        private AbsoluteDate reference;

        /** Fitted parameters. */
        private double[] fitted;

        /** Observed points. */
        private List<WeightedObservedPoint> observedPoints;

        /** Simple constructor.
         * @param secularDegree degree of polynomial secular part
         * @param pulsations pulsations of harmonic part
         */
        public  SecularAndHarmonic(final int secularDegree, final double pulsations) {
            this.secularDegree  = secularDegree;
            this.pulsations     = pulsations.clone();
            this.observedPoints = new ArrayList<WeightedObservedPoint>();
        }

        /** Reset fitting.
         * @param date reference date
         * @param initialGuess initial guess for the parameters
         * @see #getReferenceDate()
         */
        public  void resetFitting(final AbsoluteDate date, final double initialGuess) {
            reference = date;
            fitted    = initialGuess.clone();
            observedPoints.clear();
        }

        /** Add a fitting point.
         * @param date date of the point
         * @param osculatingValue osculating value
         */
        public  void addPoint(final AbsoluteDate date, final double osculatingValue) {
            observedPoints.add(new WeightedObservedPoint(1.0, date.durationFrom(reference), osculatingValue));
        }

        /** Get the reference date.
         * @return reference date
         * @see #resetFitting(AbsoluteDate, double...)
         */
        public  AbsoluteDate getReferenceDate() {
            return reference;
        }

        /** Get an upper bound of the fitted harmonic amplitude.
         * @return upper bound of the fitted harmonic amplitude
         */
        public  double getHarmonicAmplitude() {
            double amplitude = 0;
            for (int i = 0; i < pulsations.length; ++i) {
                amplitude += FastMath.hypot(fitted[secularDegree + 2 * i + 1],
                                            fitted[secularDegree + 2 * i + 2]);
            }
            return amplitude;
        }

        /** Fit parameters.
         * @see #getFittedParameters()
         */
        public  void fit() {

            final AbstractCurveFitter fitter = new AbstractCurveFitter() {
                /** {@inheritDoc} */
                @Override
                protected  LeastSquaresProblem getProblem(final Collection<WeightedObservedPoint> observations) {
                    // Prepare least-squares problem.
                    final int len = observations.size();
                    final double[] target  = new double[len];
                    final double[] weights = new double[len];

                    int i = 0;
                    for (final WeightedObservedPoint obs : observations) {
                        target[i]  = obs.getY();
                        weights[i] = obs.getWeight();
                        ++i;
                    }

                    final AbstractCurveFitter.TheoreticalValuesFunction model =
                            new AbstractCurveFitter.TheoreticalValuesFunction(new LocalParametricFunction(), observations);

                    // build a new least squares problem set up to fit a secular and harmonic curve to the observed points
                    return new LeastSquaresBuilder().
                            maxEvaluations(Integer.MAX_VALUE).
                            maxIterations(Integer.MAX_VALUE).
                            start(fitted).
                            target(target).
                            weight(new DiagonalMatrix(weights)).
                            model(model.getModelFunction(), model.getModelFunctionJacobian()).
                            build();

                }
            };

            fitted = fitter.fit(observedPoints);

        }

        /** Local parametric function used for fitting. */
        private class LocalParametricFunction implements ParametricUnivariateFunction {

            /** {@inheritDoc} */
            public  double value(final double x, final double... parameters) {
                return truncatedValue(secularDegree, pulsations.length, x, parameters);
            }

            /** {@inheritDoc} */
            public  double[] gradient(final double x, final double... parameters) {
                final double[] gradient = new double[secularDegree + 1 + 2 * pulsations.length];

                // secular part
                double xN = 1.0;
                for (int i = 0; i <= secularDegree; ++i) {
                    gradient[i] = xN;
                    xN *= x;
                }

                // harmonic part
                for (int i = 0; i < pulsations.length; ++i) {
                    gradient[secularDegree + 2 * i + 1] = FastMath.cos(pulsations[i] * x);
                    gradient[secularDegree + 2 * i + 2] = FastMath.sin(pulsations[i] * x);
                }

                return gradient;
            }

        }

        /** Get a copy of the last fitted parameters.
         * @return copy of the last fitted parameters.
         * @see #fit()
         */
        public  double[] getFittedParameters() {
            return fitted.clone();
        }

        /** Get fitted osculating value.
         * @param date current date
         * @return osculating value at current date
         */
        public  double osculatingValue(final AbsoluteDate date) {
            return truncatedValue(secularDegree, pulsations.length,
                                  date.durationFrom(reference), fitted);
        }

        /** Get fitted osculating derivative.
         * @param date current date
         * @return osculating derivative at current date
         */
        public  double osculatingDerivative(final AbsoluteDate date) {
            return truncatedDerivative(secularDegree, pulsations.length,
                                       date.durationFrom(reference), fitted);
        }

        /** Get fitted osculating second derivative.
         * @param date current date
         * @return osculating second derivative at current date
         */
        public  double osculatingSecondDerivative(final AbsoluteDate date) {
            return truncatedSecondDerivative(secularDegree, pulsations.length,
                                             date.durationFrom(reference), fitted);
        }

        /** Get mean value, truncated to first components.
         * @param date current date
         * @param degree degree of polynomial secular part to consider
         * @param harmonics number of harmonics terms to consider
         * @return mean value at current date
         */
        public  double meanValue(final AbsoluteDate date, final int degree, final int harmonics) {
            return truncatedValue(degree, harmonics, date.durationFrom(reference), fitted);
        }

        /** Get mean derivative, truncated to first components.
         * @param date current date
         * @param degree degree of polynomial secular part to consider
         * @param harmonics number of harmonics terms to consider
         * @return mean derivative at current date
         */
        public  double meanDerivative(final AbsoluteDate date, final int degree, final int harmonics) {
            return truncatedDerivative(degree, harmonics, date.durationFrom(reference), fitted);
        }

        /** Approximate an already fitted model to polynomial only terms.
         * <p>
         * This method is mainly used in order to combine the large amplitude long
         * periods with the secular part as a new approximate polynomial model over
         * some time range. This should be used rather than simply extracting the
         * polynomial coefficients from {@link #getFittedParameters()} when some
         * periodic terms amplitudes are large (for example Sun resonance effects
         * on local solar time in sun synchronous orbits). In theses cases, the pure
         * polynomial secular part in the coefficients may be far from the mean model.
         * </p>
         * @param combinedDegree desired degree for the combined polynomial
         * @param combinedReference desired reference date for the combined polynomial
         * @param meanDegree degree of polynomial secular part to consider
         * @param meanHarmonics number of harmonics terms to consider
         * @param start start date of the approximation time range
         * @param end end date of the approximation time range
         * @param step sampling step
         * @return coefficients of the approximate polynomial (in increasing degree order),
         * using the user provided reference date
         */
        public  double[] approximateAsPolynomialOnly(final int combinedDegree, final AbsoluteDate combinedReference,
                                                    final int meanDegree, final int meanHarmonics,
                                                    final AbsoluteDate start, final AbsoluteDate end,
                                                    final double step) {
            final List<WeightedObservedPoint> points = new ArrayList<WeightedObservedPoint>();
            for (AbsoluteDate date = start; date.compareTo(end) < 0; date = date.shiftedBy(step)) {
                points.add(new WeightedObservedPoint(1.0,
                                                     date.durationFrom(combinedReference),
                                                     meanValue(date, meanDegree, meanHarmonics)));
            }
            return PolynomialCurveFitter.create(combinedDegree).fit(points);
        }

        /** Get mean second derivative, truncated to first components.
         * @param date current date
         * @param degree degree of polynomial secular part
         * @param harmonics number of harmonics terms to consider
         * @return mean second derivative at current date
         */
        public  double meanSecondDerivative(final AbsoluteDate date, final int degree, final int harmonics) {
            return truncatedSecondDerivative(degree, harmonics, date.durationFrom(reference), fitted);
        }

        /** Get value truncated to first components.
         * @param degree degree of polynomial secular part
         * @param harmonics number of harmonics terms to consider
         * @param time time parameter
         * @param parameters models parameters (must include all parameters,
         * including the ones ignored due to model truncation)
         * @return truncated value
         */
        private  double truncatedValue(final int degree, final int harmonics,
                                      final double time, final double... parameters) {

            double value = 0;

            // secular part
            double tN = 1.0;
            for (int i = 0; i <= degree; ++i) {
                value += parameters[i] * tN;
                tN    *= time;
            }

            // harmonic part
            for (int i = 0; i < harmonics; ++i) {
                value += parameters[secularDegree + 2 * i + 1] * FastMath.cos(pulsations[i] * time) +
                         parameters[secularDegree + 2 * i + 2] * FastMath.sin(pulsations[i] * time);
            }

            return value;

        }

        /** Get derivative truncated to first components.
         * @param degree degree of polynomial secular part
         * @param harmonics number of harmonics terms to consider
         * @param time time parameter
         * @param parameters models parameters (must include all parameters,
         * including the ones ignored due to model truncation)
         * @return truncated derivative
         */
        private  double truncatedDerivative(final int degree, final int harmonics,
                                           final double time, final double... parameters) {

            double derivative = 0;

            // secular part
            double tN = 1.0;
            for (int i = 1; i <= degree; ++i) {
                derivative += i * parameters[i] * tN;
                tN         *= time;
            }

            // harmonic part
            for (int i = 0; i < harmonics; ++i) {
                derivative += pulsations[i] * (-parameters[secularDegree + 2 * i + 1] * FastMath.sin(pulsations[i] * time) +
                                                parameters[secularDegree + 2 * i + 2] * FastMath.cos(pulsations[i] * time));
            }

            return derivative;

        }

        /** Get second derivative truncated to first components.
         * @param degree degree of polynomial secular part
         * @param harmonics number of harmonics terms to consider
         * @param time time parameter
         * @param parameters models parameters (must include all parameters,
         * including the ones ignored due to model truncation)
         * @return truncated second derivative
         */
        private  double truncatedSecondDerivative(final int degree, final int harmonics,
                                                 final double time, final double... parameters) {

            double d2 = 0;

            // secular part
            double tN = 1.0;
            for (int i = 2; i <= degree; ++i) {
                d2 += (i - 1) * i * parameters[i] * tN;
                tN *= time;
            }

            // harmonic part
            for (int i = 0; i < harmonics; ++i) {
                d2 += -pulsations[i] * pulsations[i] *
                      (parameters[secularDegree + 2 * i + 1] * FastMath.cos(pulsations[i] * time) +
                       parameters[secularDegree + 2 * i + 2] * FastMath.sin(pulsations[i] * time));
            }

            return d2;

        }

    }

}


 
