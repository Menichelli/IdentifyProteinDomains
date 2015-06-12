/**
 * 
 */
package tools.parser;

import global.Global;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * @author christophe
 *
 */
public class RefProtParser {

	private RefProtParser() {}

	/**
	 * 
	 * @param protNames format: Name_Species
	 * @return
	 */
	public static Map<String,String> getProteinSequences(Set<String> protNames) {
		Map<String,String> ret = new HashMap<String, String>();

		Set<String> missingProts = new HashSet<String>(protNames);

		try (BufferedReader br = new BufferedReader(new FileReader(Global.REFPROT_PATH))) {
			String sCurrentLine,protName="",seq="";
			boolean needed = false;

			while(true) {
				sCurrentLine = br.readLine();
				if(sCurrentLine==null) {
					if(needed) ret.put(protName, seq);
					break;
				}
				if(sCurrentLine.startsWith(">")) {
					if(needed) ret.put(protName, seq);
					protName = sCurrentLine.split("RepID=")[1];
					needed = missingProts.contains(protName);
					seq = "";
				} else {
					if(needed) seq+=sCurrentLine;
				}
			}

		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}


		return ret;
	}

	/**
	 * Renvoie la liste permettant de convertir un ClusterID en ProtID
	 * @return
	 */
	public static Map<String,String> getClusterProtMapping() {
		Map<String,String> ret = new HashMap<String, String>();

		try (BufferedReader br = new BufferedReader(new FileReader(Global.REFPROT_PATH))) {
			String sCurrentLine,protName="",clusterName="";
			
			while(true) {
				sCurrentLine = br.readLine();
				if(sCurrentLine==null) {
					break;
				}
				if(sCurrentLine.startsWith(">")) {
					//>UniRef50_A0A009JBL3 DNA topoisomerase (ATP-hydrolyzing) n=9 RepID=A0A009JBL3_ACIBA
					clusterName = sCurrentLine.split("\\s+")[0].substring(1);
					protName = sCurrentLine.split("RepID=")[1];
					ret.put(clusterName, protName);
				} else {
					continue;
				}
			}

		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		return ret;
	}

}
