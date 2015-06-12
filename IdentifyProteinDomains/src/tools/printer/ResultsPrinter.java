/**
 * 
 */
package tools.printer;

import global.Global;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;

import model.BlastHit;
import model.PutativeDomain;

/**
 * @author christophe
 *
 */
public class ResultsPrinter {

	private static ResultsPrinter instance;
	private BufferedWriter writer;

	private ResultsPrinter() throws IOException {
		File f = new File(Global.RESULTS_PATH);
		if(!f.exists()) f.createNewFile();
		FileWriter fw = new FileWriter(f);
		writer = new BufferedWriter(fw);
		writer.write("domainValidated nbSequences bestEvalue meanEvalue domainValidating logpvalue validatedByPfam\n");
	}

	public static ResultsPrinter getInstance() throws Exception {
		instance = (instance==null)? new ResultsPrinter():instance;
		return instance;
	}

	public void close() throws Exception {
		writer.close();
		instance = null;
	}

	public void addEntry(PutativeDomain domain, String domainValidating, double score, boolean validatedByPfam) throws Exception {
		String domainValidated = domain.getIdentifier();
		Set<BlastHit> bhs = domain.getHitsCoveringResidue(domain.getBestPosition());
		int nbSequences = bhs.size();
		double bestEvalue = getBestEvalue(bhs);
		double meanEvalue = getMeanEvalue(bhs);

		String msg = domainValidated+" "+nbSequences+" "+bestEvalue+" "+meanEvalue+" "+domainValidating+" "+score+" "+validatedByPfam;
		writer.append(msg+"\n");
		writer.flush();
	}

	
	
	//TOOLS
	
	private static double getBestEvalue(Set<BlastHit> set) {
		double ret=0;

		for(BlastHit bh : set) {
			if(ret==0) ret = bh.geteValue();
			else {
				if(ret>bh.geteValue()) ret = bh.geteValue();
			}
		}

		return ret;
	}

	private static double getMeanEvalue(Set<BlastHit> set) {
		double ret = 0;

		for(BlastHit bh : set) {
			ret+=(bh.geteValue()/set.size());
		}

		return ret;
	}

}
