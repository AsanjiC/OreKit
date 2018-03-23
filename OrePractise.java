package org.orekit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Locale;
import java.util.Scanner;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.orekit.bodies.FieldGeodeticPoint;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.forces.ForceModel;
import org.orekit.forces.drag.DragForce;
import org.orekit.forces.drag.DragSensitive;
import org.orekit.forces.drag.IsotropicDrag;
import org.orekit.forces.drag.atmosphere.Atmosphere;
import org.orekit.forces.drag.atmosphere.SimpleExponentialAtmosphere;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.propagation.AbstractPropagator;
import org.orekit.propagation.analytical.AbstractAnalyticalPropagator;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.Propagator;
public class OrePractise {
	 
	public static void main(String[] args) throws IOException {
        try {

        		//Set HDD a main Directory
            File home       = new File(System.getProperty("user.home"));
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
            DataProvidersManager manager = DataProvidersManager.getInstance();
            manager.addProvider(new DirectoryCrawler(orekitData));
            /*File file = new File("ISSTLE.txt");
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            //StringBuffer stringBuffer = new StringBuffer();*/
            String line1 = "1 25544U 98067A   18080.51488668  .00001810  00000-0  34564-4 0  9999";
            String line2 = "2 25544  51.6407  95.7783 0001826 226.4922 252.5706 15.54122438104881";
            
            //fr.close();
            
            // create pseudo-orbit
            TLE tle = new TLE(line1, line2);
            // create propagator
            TLEPropagator propagator = TLEPropagator.selectExtrapolator(tle);
            
            // Inertial frame 
            Frame inertialFrame = propagator.getFrame();

            // Initial date in UTC time scale
            TimeScale utc = TimeScalesFactory.getUTC();
            AbsoluteDate initialDate = new AbsoluteDate(2017, 12, 20, 23, 30, 00.000, utc);

            // gravitation coefficient
            /*double mu =  3.986004415e+14;
            
            // Initial orbit parameters
            double a = 6600000; // semi major axis in meters
            double e = 0.0167; // eccentricity
            double i = FastMath.toRadians(7); // inclination
            double omega = FastMath.toRadians(180); // perigee argument
            double raan = FastMath.toRadians(261); // right ascension of ascending node
            double lM = 0; // mean anomaly
            */	
            // Orbit construction as Keplerian
            Orbit initialOrbit = propagator.propagateOrbit(initialDate);
            
            //initial StatenDefinition
            SpacecraftState initialState = new SpacecraftState(initialOrbit);
            
            // Step integrator 
            final double minStep = 0.001;
            final double maxstep = 2000.0;
            final double positionTolerance = 10.0;
            final OrbitType propagationType = OrbitType.KEPLERIAN;
            final double[][] tolerances = NumericalPropagator.tolerances(positionTolerance, initialOrbit, propagationType);
            AdaptiveStepsizeIntegrator integrator = new DormandPrince853Integrator(minStep, maxstep, tolerances[0], tolerances[1]);

            // Propagator
            //NumericalPropagator propagator = new NumericalPropagator(integrator);
            
            //propagator.setOrbitType(propagationType);

            // Force Model
            final double ae = Constants.GRS80_EARTH_EQUATORIAL_RADIUS;
            OneAxisEllipsoid earth = new OneAxisEllipsoid(ae, Constants.WGS84_EARTH_FLATTENING, FramesFactory.getEME2000());
            earth.setAngularThreshold(1.e-6);
            Atmosphere atm = new SimpleExponentialAtmosphere(earth, 4.e-13, 200000.0, 8500.0);
            final double cd = 2.0;
            final double sf = 5.0;
            ForceModel dragNUM = new DragForce(atm, new IsotropicDrag(sf, cd));

            // Add force model to the propagator
            //propagator.addForceModel(dragNUM);

            // Set up initial state in the propagator
            //propagator.setInitialState(initialState);

            // Set up operating mode for the propagator as master mode
            // with fixed step(seconds) and specialized step handler
            propagator.setMasterMode(10., new stepHandler());

            // Extrapolate from the initial to the final date
            //SpacecraftState finalState = propagator.propagate(initialDate.shiftedBy(10.));
            //Updated to Propagate for 1 Month
            SpacecraftState finalState = propagator.propagate(initialDate.shiftedBy(7.884e+6));
            
            KeplerianOrbit o = (KeplerianOrbit) OrbitType.KEPLERIAN.convertType(finalState.getOrbit());
            System.out.print(o.initPVCoordinates().getPosition()+"\n");
            
            //These NextTwo Lines is How I Convert the coordinates to get Long. Lat & Altitude 
            FieldGeodeticPoint<DerivativeStructure> gp = earth.transform(o.initPVCoordinates(),inertialFrame,initialDate.shiftedBy(7.884e+6));
            System.out.print(gp.toString() +"\n");
            System.out.println(propagator.getTLE());
            System.out.format(Locale.US, "Final state:%n%s %12.3f %10.8f %10.6f %10.6f %10.6f %10.6f%n",
                              finalState.getDate(),
                              o.getA(), o.getE(),
                              FastMath.toDegrees(o.getI()),
                              FastMath.toDegrees(o.getPerigeeArgument()),
                              FastMath.toDegrees(o.getRightAscensionOfAscendingNode()),
                              FastMath.toDegrees(o.getTrueAnomaly()));  
        } catch (OrekitException oe) {
            System.err.println(oe.getMessage());
        }
    }
    
	
    
    // Step Handler Class used to print on the output stream at the given step.
    private static class stepHandler implements OrekitFixedStepHandler {
    	
    	
    	String line1 = "1 25544U 98067A   18080.51488668  .00001810  00000-0  34564-4 0  9999";
        String line2 = "2 25544  51.6407  95.7783 0001826 226.4922 252.5706 15.54122438104881";
    	TLEPropagator propagator;
    	
    	
    	final double ae = Constants.GRS80_EARTH_EQUATORIAL_RADIUS;
        OneAxisEllipsoid earth = new OneAxisEllipsoid(ae, Constants.WGS84_EARTH_FLATTENING, FramesFactory.getEME2000());
    	

        private stepHandler(){

        }
        /*
        public  void init(final SpacecraftState s0, final AbsoluteDate t, final double step) {
            System.out.println("          date                a           e" +
                               "           i         \u03c9          \u03a9" +
                               "          \u03bd");
        }
        */
        //Method to create a propagator from a Two Line Element
        public TLEPropagator createPropagator(final String line1, final String line2) throws OrekitException 
        {
                // create pseudo-orbit
                TLE tle = new TLE(line1, line2);
                // create propagator
                TLEPropagator propagator = TLEPropagator.selectExtrapolator(tle);
                return propagator;
        }
        
        public  void handleStep(SpacecraftState currentState, boolean isLast) throws OrekitException {
        	propagator = createPropagator(line1,line2);
        	Frame inertialFrame = propagator.getFrame();
        	TimeScale utc = TimeScalesFactory.getUTC();
        	AbsoluteDate initialDate = new AbsoluteDate(2017, 12, 20, 23, 30, 00.000, utc);
            KeplerianOrbit o = (KeplerianOrbit) OrbitType.KEPLERIAN.convertType(currentState.getOrbit());
            //System.out.print(o.initPVCoordinates().getPosition()+"\n"); 
            System.out.println(propagator.getPVCoordinates(initialDate.shiftedBy(7.884e+6)));
            FieldGeodeticPoint<DerivativeStructure> gp = earth.transform(o.initPVCoordinates(),inertialFrame,initialDate.shiftedBy(7.884e+6));
            System.out.print(o.getDate() + " " + gp.toString() +"\n");
            
            
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
  
    }
    


 
