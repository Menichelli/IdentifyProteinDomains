/**
 * 
 */
package module;

import global.Global;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import model.BlastHit;
import model.CouplePfamPutative;
import model.PfamFamily;
import model.PutativeDomain;
import model.ValidatedDomain;
import tools.Collection;
import tools.CoocScoringModule;
import tools.CoupleGenerator;
import tools.HitsGathering;
import tools.PutativeShuffler;
import tools.parser.BlastResultsParser;
import tools.parser.PfamParser;
import tools.printer.ConservationStatsPrinter;
import tools.printer.FastaPrinter;
import tools.printer.StatsPrinter;

/**
 * @author christophe
 *
 */
public class PfamValidationModule extends AbstractValidationModule {

	private final static double hitsEvalueMax = 1e-2;

	@Override
	public void run() {
		try {
			//Step 1: Recuperer tous les Pfam
			if(Global.VERBOSE) System.out.println("Parsing Pfam families...");
			Set<PfamFamily> pfamFamilies = PfamParser.getFamilies();
			Map<String,PfamFamily> mapIdPfam = new HashMap<String, PfamFamily>();
			for(PfamFamily pf : pfamFamilies) {
				mapIdPfam.put(pf.getFamilyName(), pf);
			}
			if(Global.VERBOSE) System.out.println("Found "+pfamFamilies.size()+" different families.");

			//Step 2: Recuperer tous les hits avec une bonne p-valeur
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

			//Step 4: Construire tous les couples possibles
			if(Global.VERBOSE) System.out.println("Testing all couples...");
			int nbCouplesTotal = 0, nbCouplesRetained = 0;
			int nbProtTreated = 0;
			for(String protName : putativeDomainsByProt.keySet()) {
				Set<CouplePfamPutative> couples = CoupleGenerator.getCouplePfamPutative(putativeDomainsByProt.get(protName), pfamFamilies);
				nbCouplesTotal+=couples.size();
				for(CouplePfamPutative couple : couples) {
					String pfamFamilyName = couple.getPfam().getFamilyName();
					String putativeDomainIdentifier = couple.getPutative().getIdentifier();
					int nbProtPfam = couple.getPfam().getAllProteinNames().size();
					Set<String> protsCoveringThePutativeDomain = couple.getPutative().getProteinsCoveringResidue(couple.getPutative().getBestPosition());
					int nbProtPutativeDomain = protsCoveringThePutativeDomain.size();
					int nbProtIntersec = Collection.intersectionSize(couple.getPfam().getAllProteinNames(), protsCoveringThePutativeDomain);

					if(nbProtIntersec >= Global.NB_SEQ_INTERSECT) {
						nbCouplesRetained++;
						StatsPrinter.getInstance(Global.STATS_PFPT_PATH).addEntry(pfamFamilyName, putativeDomainIdentifier, nbProtPfam, nbProtPutativeDomain, nbProtIntersec);
					}
				}
				nbProtTreated++;
				if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.print("\r"+nbProtTreated+"/"+putativeDomainsByProt.keySet().size()+" proteins tested.");
			}
			if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.println();
			StatsPrinter.getInstance(Global.STATS_PFPT_PATH).close();
			if(Global.VERBOSE) System.out.println("Generated "+nbCouplesTotal+" couples Pfam-Putative. "+nbCouplesRetained+" couples with at least "+Global.NB_SEQ_INTERSECT+" proteins in common.");

			//Step 5: Tester le score de cooc et garder les meilleurs
			if(Global.VERBOSE) System.out.print("Computing coocurrence scores...");
			Set<ValidatedDomain> validatedDomains = CoocScoringModule.compute(Global.STATS_PFPT_PATH, Global.R_RESULTS_PFPT_PATH);
			if(Global.VERBOSE) System.out.println("done.");
			if(Global.KEEPONLYBESTPVALUE) {
				Map<String,ValidatedDomain> keepBestPvalue = new HashMap<String, ValidatedDomain>();
				ValidatedDomain tmpVD;
				for(ValidatedDomain vd : validatedDomains) {
					tmpVD = keepBestPvalue.get(vd.getIdentifierValidatedDomain());
					if(tmpVD!=null) {
						if(vd.getScore()<tmpVD.getScore()) {
							keepBestPvalue.put(vd.getIdentifierValidatedDomain(), vd);
						}
					} else {
						keepBestPvalue.put(vd.getIdentifierValidatedDomain(), vd);
					}
					validatedDomains.clear();
					for(String s : keepBestPvalue.keySet()) {
						validatedDomains.add(keepBestPvalue.get(s));
					}
				}
			}
			
			//Step 6: Print le modele de chaque validation
			if(Global.VERBOSE) System.out.print("Initializing the FastaPrinter...");
			Set<String> initSet = new HashSet<String>();
			for(ValidatedDomain vd : validatedDomains) {
				initSet.addAll(mapIdPfam.get(vd.getIdentifierValidatingDomain()).getAllProteinNames());
			}
			FastaPrinter.getInstance().init(initSet);
			if(Global.VERBOSE) System.out.println("Ready.");
			
			if(Global.VERBOSE) System.out.println("Printing...");
			int domPrinted = 0;
			int nbPrinted;
			for(ValidatedDomain vd : validatedDomains) {
				nbPrinted = FastaPrinter.getInstance().printFasta(mapIdPutativeDomain.get(vd.getIdentifierValidatedDomain()), mapIdPfam.get(vd.getIdentifierValidatingDomain()).getAllProteinNames());
				ConservationStatsPrinter.getInstance("PfamValidationConservation.dat").addEntry(nbPrinted, mapIdPutativeDomain.get(vd.getIdentifierValidatedDomain()).getBlastHits().size());
				domPrinted++;
				if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.print("\r> "+domPrinted);
			}
			ConservationStatsPrinter.getInstance("PfamValidationConservation.dat").close();
			if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.println();
			if(Global.VERBOSE) System.out.println("Printing done.");
			
			//Step 7: Estimer le fdr
			if(Global.VERBOSE) System.out.println("Computing FDR...");
			estimateFDR(putativeDomainsByProt, pfamFamilies, validatedDomains.size());
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private double estimateFDR(final Map<String,Set<PutativeDomain>> putativeDomainsByProt, final Set<PfamFamily> pfamFamilies, int nbCertificationObtained) throws Exception {
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
					int nbProtIntersec = Collection.intersectionSize(couple.getPfam().getAllProteinNames(), protsCoveringThePutativeDomain);

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
				Set<ValidatedDomain> validatedDomains = CoocScoringModule.compute(Global.FDR_TMP_PATH+"1", Global.FDR_TMP_PATH+"2");
				Set<String> vDomains = new HashSet<String>();
				for(ValidatedDomain v : validatedDomains) {
					vDomains.add(v.getIdentifierValidatedDomain());
				}
				//currentCertification = vDomains.size(); //compte 1 pour chaque domaine
				currentCertification = validatedDomains.size(); //compte 1 pour chaque validations
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
		
		double fdr = totalCertification/Global.FDR_NB_REPEATS/nbCertificationObtained;
		if(Global.VERBOSE) {
			String s = String.format("%d min, %d sec",
					TimeUnit.MILLISECONDS.toMinutes(totalTimeElapsed),
					TimeUnit.MILLISECONDS.toSeconds(totalTimeElapsed) -
					TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(totalTimeElapsed))
					);
			System.out.println("progress: 100%, finished in "+s+".                             ");
			System.out.println("FDREstimator found a total of " + totalCertification + " certifications over " + Global.FDR_NB_REPEATS + " repeats.");
			System.out.println("Process found "+nbCertificationObtained+" certifications.");
			System.out.println("FDR is estimated at: "+fdr);
		}

		return fdr;
	}



}
