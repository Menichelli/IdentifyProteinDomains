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

import model.BlastHit;
import model.CouplePfamPutative;
import model.PfamFamily;
import model.PutativeDomain;
import model.ValidatedDomain;
import tools.Collection;
import tools.CoocScoringModule;
import tools.CoupleGenerator;
import tools.HitsGathering;
import tools.parser.BlastResultsParser;
import tools.parser.PfamParser;
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
				System.out.println("Found "+nbhit+" valid hits on "+hitsByProt.keySet().size()+" different proteins.");
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
				initSet.add(mapIdPutativeDomain.get(vd.getIdentifierValidatedDomain()).getQueryName()+"_"+mapIdPutativeDomain.get(vd.getIdentifierValidatedDomain()).getQuerySpecies());
			}
			FastaPrinter.getInstance().init(initSet);
			if(Global.VERBOSE) System.out.println("Ready.");
			
			for(ValidatedDomain vd : validatedDomains) {
				FastaPrinter.getInstance().printFasta(mapIdPutativeDomain.get(vd.getIdentifierValidatedDomain()), mapIdPfam.get(vd.getIdentifierValidatingDomain()).getAllProteinNames());
			}
			
			//Step 7: Estimer le fdr
			if(Global.VERBOSE) System.out.println("FDR estimator NYI");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}



}
