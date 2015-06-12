/**
 * 
 */
package tools.parser;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashSet;
import java.util.Set;

import model.ValidatedDomain;

/**
 * @author christophe
 *
 */
public class RParser {
	
	private RParser() {}
	
	public static Set<ValidatedDomain> getValidatedDomains(String path) throws Exception {
		Set<ValidatedDomain> ret = new HashSet<ValidatedDomain>();
		
		BufferedReader br = new BufferedReader(new FileReader(new File(path)));
		String sCurrentLine;
		String[] tCurrentLine;
		
//		pfamFamilyName putativeDomainIdentifier nbProtPfam nbProtPutativeDomain nbProtIntersec log(pValue)
//		CDC48_N O77313_PLAF7/427-562   1231 178  43 -193.501389248843
//		AAA O77313_PLAF7/427-562  47246 178  93 -159.792236886434
		br.readLine(); //skip the header
		double pvalue;
		while((sCurrentLine = br.readLine())!=null) {
			tCurrentLine = sCurrentLine.split("\\s+");
			if(tCurrentLine==null || sCurrentLine.equals("NULL")) break;
			if(tCurrentLine[5].equals("-Inf")) pvalue = Double.MIN_VALUE;
			else pvalue = Double.parseDouble(tCurrentLine[5]);
			ret.add(new ValidatedDomain(tCurrentLine[1], tCurrentLine[0], pvalue));
		}
		br.close();
		
		return ret;
	}
	
}
