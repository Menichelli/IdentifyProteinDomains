import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import tools.HitsGathering;
import tools.parser.BlastResultsParser;
import tools.parser.ProteomeParser;
import model.BlastHit;
import model.PutativeDomain;
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

	private final static double hitsEvalueMax = 1e-2;

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

		//Step : Parser les resultats BLAST
		if(Global.VERBOSE) System.out.println("Parsing Blast results...");
		Map<String,List<BlastHit>> hitsByProt = BlastResultsParser.getHitsByProt(hitsEvalueMax);
		if(Global.VERBOSE) {
			int nbhit = 0;
			for(String s : hitsByProt.keySet()) {
				nbhit += hitsByProt.get(s).size();
			}
			System.out.println("Found "+nbhit+" valid hits on "+hitsByProt.keySet().size()+" different proteins (max evalue: "+hitsEvalueMax+").");
		}

		//Step 3: Construire tous les domaines potentiels
		if(Global.VERBOSE) System.out.println("Gathering putative domains...");
		Map<String,Set<PutativeDomain>> putativeDomainsByProt = new HashMap<String, Set<PutativeDomain>>();
		Set<PutativeDomain> tmpset,tmp;
		for(String protName : hitsByProt.keySet()) {
			List<BlastHit> hitsOnThisProt = hitsByProt.get(protName);
			while(true) { //break si le gatherHits ne trouve pas de putative
				tmpset = HitsGathering.gatherHits(hitsOnThisProt,Global.PROTEOME_AIMED.getProteinByName(protName),putativeDomainsByProt.get(protName));
				if(!tmpset.isEmpty()) {
					tmp = putativeDomainsByProt.get(protName);
					if(tmp==null) tmp = new HashSet<PutativeDomain>();
					tmp.addAll(tmpset);
					putativeDomainsByProt.put(protName,tmp);
					for(PutativeDomain pdom : tmpset) { //sur chaque domaine trouve
						for(BlastHit h : pdom.getBlastHits()) { //sur chacun de ses hits
							hitsOnThisProt.remove(h); //pour ne pas reutiliser le meme hit deux fois
						}
					}
				} else break;
			}
		}
		if(Global.VERBOSE) {
			int nbhit = 0;
			for(String s : putativeDomainsByProt.keySet()) {
				nbhit += putativeDomainsByProt.get(s).size();
			}
			System.out.println("Found "+nbhit+" putative domains on "+putativeDomainsByProt.keySet().size()+" different proteins.");
		}
		Map<String,PutativeDomain> mapIdPutativeDomain = new HashMap<String, PutativeDomain>();
		for(String s : putativeDomainsByProt.keySet()) {
			for(PutativeDomain pd : putativeDomainsByProt.get(s)) {
				mapIdPutativeDomain.put(pd.getIdentifier(), pd);
			}
		}

		//----------------------------------------------------------------------------------//
		//----------------------------------------------------------------------------------//
		//----------------------------------------------------------------------------------//


		//Step 2: Run PfamValidation
		if(Global.VERBOSE) System.out.println("***Starting PfamValidationModule...");
		if(Global.VERBOSE) printState(putativeDomainsByProt);
		AbstractValidationModule m1 = new PfamValidationModule(putativeDomainsByProt,mapIdPutativeDomain);
		m1.start();
		m1.join();
		if(Global.VERBOSE) printState(putativeDomainsByProt);
		if(Global.VERBOSE) System.out.println("***PfamValidationModule done.\n");

		//Step 3: Run CrossValidation
		if(Global.VERBOSE) System.out.println("***Starting CrossValidationModule...");
		if(Global.VERBOSE) printState(putativeDomainsByProt);
		AbstractValidationModule m2 = new CrossValidationModule(putativeDomainsByProt,mapIdPutativeDomain);
		m2.start();
		m2.join();
		if(Global.VERBOSE) printState(putativeDomainsByProt);
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

	private static void printState(Map<String,Set<PutativeDomain>> putativeDomainsByProt) {
		int nbCertified = 0;
		int nbTotal = 0;
		for(String s : putativeDomainsByProt.keySet()) {
			for(PutativeDomain pd : putativeDomainsByProt.get(s)) {
				if(pd.hasBeenCertified()) nbCertified++;
				nbTotal++;
			}
		}
		System.out.println("** Putative domains certified: "+nbCertified+"/"+nbTotal+" **");
	}

}
