/**
 * 
 */
package tools.parser;

import global.Global;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import model.Protein;
import model.Proteome;

/**
 * @author christophe
 *
 */
public class ProteomeParser {

	private ProteomeParser() {}

	public static Proteome getProteome() {
		Proteome ret = new Proteome();
		try (BufferedReader br = new BufferedReader(new FileReader(Global.PROTEOME_PATH))) {
			String sCurrentLine,currentProtName="",currentSeq="",currentProtID="",currentSpecies="";

			while((sCurrentLine = br.readLine())!=null) {
				if(sCurrentLine.startsWith(">")) {
					if(!currentSeq.isEmpty()) {
						ret.addProtein(new Protein(currentProtName, currentProtID, currentSpecies, currentSeq));
					}
					currentProtName = "";
					currentSeq = "";
					currentProtID = "";
					currentSpecies = "";
					currentProtName = sCurrentLine.split("\\s+")[0].split("\\|")[2].split("_")[0];
					currentProtID = sCurrentLine.split("\\s+")[0].split("\\|")[1];
					currentSpecies = sCurrentLine.split("\\s+")[0].split("\\|")[2].split("_")[1];
				} else {
					currentSeq += sCurrentLine;
				}
			}
			ret.addProtein(new Protein(currentProtName, currentProtID, currentSpecies, currentSeq));
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return ret;
	}
	
	public static Map<String, Integer> getProteinsSize() {
		Map<String, Integer> ret = new HashMap<String, Integer>();
		try (BufferedReader br = new BufferedReader(new FileReader(Global.PROTEOME_PATH))) {
			String sCurrentLine,currentProtName="",currentSeq="";

			while((sCurrentLine = br.readLine())!=null) {
				if(sCurrentLine.startsWith(">")) {
					if(!currentSeq.isEmpty()) {
						ret.put(currentProtName, currentSeq.length());
					}
					currentProtName = "";
					currentSeq = "";
					currentProtName = sCurrentLine.split("\\s+")[0].split("\\|")[2].split("_")[0];
				} else {
					currentSeq += sCurrentLine;
				}
			}
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return ret;
	}
	
}
