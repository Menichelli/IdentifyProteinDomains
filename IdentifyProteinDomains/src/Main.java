import java.util.Set;
import java.util.concurrent.TimeUnit;

import tools.parser.PfamParser;
import tools.parser.ProteomeParser;
import model.PfamFamily;
import module.PfamValidationModule;
import global.Global;

/**
 * 
 */

/**
 * @author christophe
 */
public class Main {
	
	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();
		
		//Step 1: Initialize properties
		if(args.length==1) {
			Global.PATHPROPERTIES = args[0];
		} else Global.PATHPROPERTIES = "ppblast.properties";
		System.out.print("Initialization...");
		Global.init();
		Global.PROTEOME_AIMED = ProteomeParser.getProteome();
		System.out.println(" done.");
		
		if(Global.VERBOSE) System.out.println("Found "+Global.PROTEOME_AIMED.getProteins().size()+" proteins in the proteome.");
		
		//Step 2: Run PfamValidation
		if(Global.VERBOSE) System.out.println("Starting PfamValidationModule...");
		new PfamValidationModule().start();
		if(Global.VERBOSE) System.out.println("PfamValidationModule done.");
		
		//Step 3: Run CrossValidation
		//TODO
		
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
