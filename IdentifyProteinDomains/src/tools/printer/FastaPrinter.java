/**
 * 
 */
package tools.printer;

import global.Global;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import tools.parser.RefProtParser;
import model.BlastHit;
import model.Protein;
import model.PutativeDomain;

/**
 * @author christophe
 *
 */
public class FastaPrinter {

	private static FastaPrinter instance;
	private static Map<String,String> proteins;

	private final static int CHUNK_SIZE = 60;

	private FastaPrinter() {
		proteins = new HashMap<String, String>();
	}

	public static FastaPrinter getInstance() {
		instance = (instance==null)?new FastaPrinter():instance;
		return instance;
	}

	public static void close() {
		instance = null;
	}

	public void init(Set<String> proteins) {
		FastaPrinter.proteins.putAll(RefProtParser.getProteinSequences(proteins));
	}

	public int printFasta(PutativeDomain domain, Set<String> allowedProteins) throws Exception {
		int ret = 0;
		
		Set<BlastHit> hits = domain.getHitsCoveringResidue(domain.getBestPosition()); //all blast hits on this domain
		Set<BlastHit> coocshits = new HashSet<BlastHit>(); //only hits found with coocs

		for(BlastHit hit : hits) {
			if(allowedProteins.contains(hit.getSubjectName()+"_"+hit.getSubjectSpecies())) {
				coocshits.add(hit);
			}
		}

		Set<String> missingProts = extractProtNames(coocshits);
		missingProts.removeAll(proteins.keySet());
		if(!missingProts.isEmpty()) {
			proteins.putAll(RefProtParser.getProteinSequences(missingProts));
		}

		File f = null;
		boolean created = false;
		int index = 0;
		while(!created) {
			f = new File(Global.FASTA_DIR+"/"+domain.getFileIdentifier()+"."+index+++".fasta");
			if(!f.exists()) {
				f.createNewFile();
				created = true;
			}
		}

		FileWriter fw = new FileWriter(f);
		BufferedWriter writer = new BufferedWriter(fw);
		writer.write("");

		//print the plasmodium sequence first
		writer.append(">"+domain.getQueryName()+"_"+domain.getQuerySpecies()+"/"+domain.getDomainStart()+"-"+domain.getDomainEnd()+"\n");
		Protein protPlaf7 = Global.PROTEOME_AIMED.getProteinByName(domain.getQueryName()+"_"+domain.getQuerySpecies());
		String seq = protPlaf7.getSequence().substring(domain.getDomainStart()-1, domain.getDomainEnd()-1);
		for(String s : splitByNumber(seq, CHUNK_SIZE)) {
			writer.append(s+"\n");
		}
		ret++;
		
		//print hits sequences
		for(BlastHit bh : coocshits) {
			writer.append(">"+bh.getSubjectName()+"_"+bh.getSubjectSpecies()+"/"+bh.getsStart()+"-"+bh.getsEnd()+"\n");
			seq = proteins.get(bh.getSubjectName()+"_"+bh.getSubjectSpecies()).substring(bh.getsStart()-1, bh.getsEnd()-1);
			for(String s : splitByNumber(seq, CHUNK_SIZE)) {
				writer.append(s+"\n");
			}
			ret++;
		}

		writer.close();
		return ret;
	}

	/**
	 * Permet d'extraire le nom des proteines concernees par ces hits
	 * @param hits
	 * @return
	 */
	private Set<String> extractProtNames(Set<BlastHit> hits) {
		Set<String> ret = new HashSet<String>();
		for(BlastHit bh : hits) {
			ret.add(bh.getSubjectName()+"_"+bh.getSubjectSpecies());
		}
		return ret;
	}

	/**
	 * Split une string suivant une taille de chunk
	 * @param s
	 * @param chunkSize
	 * @return
	 */
	private static String[] splitByNumber(String s, int chunkSize){
		int chunkCount = (s.length() / chunkSize) + (s.length()%chunkSize == 0 ? 0 : 1);
		String[] returnVal = new String[chunkCount];
		for(int i=0;i<chunkCount;i++){
			returnVal[i] = s.substring(i*chunkSize, Math.min((i+1)*chunkSize, s.length()));
		}
		return returnVal;
	}

}
