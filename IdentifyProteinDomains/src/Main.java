import java.util.concurrent.TimeUnit;

import tools.parser.ProteomeParser;
import module.AbstractValidationModule;
import module.CrossValidationModule;
import module.PfamValidationModule;
import global.Global;

/**
 * 
 */

/**
 * @author christophe
 */
public class Main {

	public static void main(String[] args) throws Exception {
		long startTime = System.currentTimeMillis();

		//Step 1: Initialize properties
		if(args.length==1) {
			Global.PATHPROPERTIES = args[0];
		} else Global.PATHPROPERTIES = "ppblast.properties";
		if(Global.VERBOSE) System.out.println("Initialization...");
		Global.init();
		Global.PROTEOME_AIMED = ProteomeParser.getProteome();
		if(Global.VERBOSE) System.out.println("Init done.");
		if(Global.VERBOSE) System.out.println("Found "+Global.PROTEOME_AIMED.getProteins().size()+" proteins in the proteome.");
		if(Global.VERBOSE) System.out.println();

		//Step 2: Run PfamValidation
		if(Global.VERBOSE) System.out.println("***Starting PfamValidationModule...");
		AbstractValidationModule m1 = new PfamValidationModule();
		m1.start();
		m1.join();
		if(Global.VERBOSE) System.out.println("***PfamValidationModule done.\n");
		
		//Step 3: Run CrossValidation
		if(Global.VERBOSE) System.out.println("***Starting CrossValidationModule...");
		AbstractValidationModule m2 = new CrossValidationModule();
		m2.start();
		m2.join();
		if(Global.VERBOSE) System.out.println("***CrossValidationModule done.");
		
		//Final step: Print execution time
		long stopTime = System.currentTimeMillis();
		long millis = stopTime-startTime;
		String s = String.format("%d min, %d sec",
				TimeUnit.MILLISECONDS.toMinutes(millis),
				TimeUnit.MILLISECONDS.toSeconds(millis) -
				TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
				);
		System.out.println("\nExecution Time: "+s);
	}

}
