/**
 * 
 */
package tools;

import global.Global;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import tools.parser.RParser;
import tools.printer.StatsPrinter;
import model.CouplePfamPutative;
import model.PfamFamily;
import model.PutativeDomain;
import model.ValidatedDomain;

import com.google.common.collect.Sets;

/**
 * @author christophe
 *
 */
public class FDREstimator {

	private FDREstimator() {}

	public static double estimateFDR(final Map<String,Set<PutativeDomain>> putativeDomainsByProt, final Set<PfamFamily> pfamFamilies, int nbPFPTobtained) throws Exception {
		double totalCertification = 0;
		
		int max=Global.FDR_NB_REPEATS,percent=0;
		
		long totalTimeElapsed = 0;
		long currentStartTime,meanTimeByRun,timeRemaining;
		
		if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.print("progress: "+percent+"%, time remaining about: estimating...\r");
		
		for(int nbRepeats = 0; nbRepeats < max; nbRepeats++) {
			currentStartTime = System.currentTimeMillis();
			
			//Step 1: Shuffle les putative domains
			Map<String,Set<PutativeDomain>> putativeDomainsByProtShuffled = PutativeShuffler.getInstance().shuffle(putativeDomainsByProt);

			//Step 2: Genere tous les couples Pfam-PutativeDomain et envoie les info au StatsPrinter
			boolean atLeastOneEntry = false;
			for(String protName : putativeDomainsByProtShuffled.keySet()) {
				Set<CouplePfamPutative> couples = CoupleGenerator.getCouplePfamPutative(putativeDomainsByProtShuffled.get(protName), pfamFamilies);
				for(CouplePfamPutative couple : couples) {
					String pfamFamilyName = couple.getPfam().getFamilyName();
					String putativeDomainIdentifier = couple.getPutative().getIdentifier();
					int nbProtPfam = couple.getPfam().getAllProteinNames().size();
					Set<String> protsCoveringThePutativeDomain = couple.getPutative().getProteinsCoveringResidue(couple.getPutative().getBestPosition());
					int nbProtPutativeDomain = protsCoveringThePutativeDomain.size();
					int nbProtIntersec = Sets.intersection(couple.getPfam().getAllProteinNames(), protsCoveringThePutativeDomain).size();

					if(nbProtIntersec >= Global.NB_SEQ_INTERSECT) {
						atLeastOneEntry|=true;
						StatsPrinter.getInstance(Global.FDR_TMP_PATH+"1").addEntry(pfamFamilyName, putativeDomainIdentifier, nbProtPfam, nbProtPutativeDomain, nbProtIntersec);
					}
				}
			}
			StatsPrinter.getInstance(Global.FDR_TMP_PATH+"1").close();
			
			int currentCertification = 0;
			if(atLeastOneEntry) {
				//Step 3: run R
				Process child1 = Runtime.getRuntime().exec("Rscript "+Global.R_SCRIPT_PATH+" "+Global.FDR_TMP_PATH+"1"+" "+Global.FDR_TMP_PATH+"2");
				int code1 = child1.waitFor();
				switch (code1) {
				case 0:
					break;
				case 1:
					throw new Exception("R script failed.");
				}

				//Step 4: gather R results
				Set<ValidatedDomain> validatedDomains = RParser.getValidatedDomains(Global.FDR_TMP_PATH+"2");
				Set<String> vDomains = new HashSet<String>();
				for(ValidatedDomain v : validatedDomains) {
					vDomains.add(v.getIdentifierValidatedDomain());
				}
				currentCertification = vDomains.size();
			}
			totalCertification += currentCertification;
			
			totalTimeElapsed += System.currentTimeMillis() - currentStartTime;
			meanTimeByRun = totalTimeElapsed / (nbRepeats+1);
			
			timeRemaining = (max-(nbRepeats+1)) * meanTimeByRun;
			
			percent = (int)((double)((double)(nbRepeats+1)/(double)max)*100);
			String s = String.format("%d min, %d sec",
					TimeUnit.MILLISECONDS.toMinutes(timeRemaining),
					TimeUnit.MILLISECONDS.toSeconds(timeRemaining) -
					TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(timeRemaining))
					);
			if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.print("progress: "+percent+"%, time remaining about: "+s+"                           \r");
		}
		
		
		
		double fdr = totalCertification/Global.FDR_NB_REPEATS/nbPFPTobtained;
		if(Global.VERBOSE) {
			String s = String.format("%d min, %d sec",
					TimeUnit.MILLISECONDS.toMinutes(totalTimeElapsed),
					TimeUnit.MILLISECONDS.toSeconds(totalTimeElapsed) -
					TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(totalTimeElapsed))
					);
			System.out.println("progress: 100%, finished in "+s+".                             ");
			System.out.println("FDREstimator found a total of " + totalCertification + " validations over " + Global.FDR_NB_REPEATS + " repeats.");
			System.out.println("Process found "+nbPFPTobtained+" validated domains.");
			System.out.println("FDR is estimated at: "+fdr);
		}

		return fdr;
	}

}
