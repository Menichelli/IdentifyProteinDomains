package tools.parser;

import global.Global;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
/**
 * 
 */


import java.util.Map;

import model.BlastHit;

/**
 * @author christophe
 *
 */
public class BlastResultsParser {

	private BlastResultsParser() {}

	public static Map<String,List<BlastHit>> getHitsByProt() throws Exception {
		return getHitsByProt(0);
	}

	public static Map<String,List<BlastHit>> getHitsByProt(double evalue_max) throws Exception {
		Map<String,List<BlastHit>> ret = new HashMap<String, List<BlastHit>>();
		List<BlastHit> hits = getHits();

		List<BlastHit> tmp;
		for(BlastHit bh : hits) {
			if(bh.geteValue() < evalue_max) {
				tmp = ret.get(bh.getQueryName()+"_"+bh.getQuerySpecies());
				if(tmp==null) tmp = new ArrayList<BlastHit>();
				tmp.add(bh);
				ret.put(bh.getQueryName()+"_"+bh.getQuerySpecies(), tmp);
			}
		}

		return ret;
	}

	public static List<BlastHit> getHits() throws Exception {
		List<BlastHit> ret = new ArrayList<BlastHit>();
		BufferedReader br = new BufferedReader(new FileReader(new File(Global.BLAST_RESULTS_PATH)));
		String sCurrentLine;
		String[] tCurrentLine;
		String clusterID, protID;

		BlastHit tmpBH;

		while((sCurrentLine = br.readLine())!=null) {
			tCurrentLine = sCurrentLine.split("\\s+");

			clusterID = tCurrentLine[1];
			protID = Global.MAPPING_CLUSTERID.get(clusterID);

			if(!protID.contains("_")) continue; //isomorph
			tmpBH = new BlastHit(
					tCurrentLine[0].split("\\|")[2].split("_")[0].intern(),
					tCurrentLine[0].split("\\|")[2].split("_")[1],
					protID.split("_")[0].intern(),
					protID.split("_")[1],
					Double.parseDouble(tCurrentLine[2])/100,
					Integer.parseInt(tCurrentLine[3]),
					Integer.parseInt(tCurrentLine[6]),
					Integer.parseInt(tCurrentLine[7]),
					Integer.parseInt(tCurrentLine[8]),
					Integer.parseInt(tCurrentLine[9]),
					Double.parseDouble(tCurrentLine[10]));
			if(!Global.EXCLUDE_SPECIES_SET.contains(tmpBH.getSubjectSpecies()) && (tmpBH.getsEnd()-tmpBH.getsStart())>=Global.SIZE_HITS_MIN) {
				ret.add(tmpBH);
			}
		}
		br.close();

		return ret;
	}

}
