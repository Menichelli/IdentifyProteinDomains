package tools.parser;

import global.Global;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashSet;
import java.util.Set;

import model.PfamFamily;

public class PfamParser {

	private PfamParser() {}

	@SuppressWarnings("unchecked")
	public static Set<PfamFamily> getFamilies() throws Exception {
		Set<PfamFamily> ret;
		if(new File(Global.PFAM_PATH+".bin").exists()) {
			ObjectInputStream ois = new ObjectInputStream(new FileInputStream(Global.PFAM_PATH+".bin"));
			ret = (Set<PfamFamily>) ois.readObject();
			ois.close();
		} else {
			ret = new HashSet<PfamFamily>();

			BufferedReader br = new BufferedReader(new FileReader(new File(Global.PFAM_PATH)));
			String sCurrentLine;

			PfamFamily currentFamily = null;

			while((sCurrentLine = br.readLine())!=null) {
				if(sCurrentLine.startsWith("#=GF ID")) { //new family
					if(currentFamily!=null) ret.add(currentFamily);
					currentFamily = new PfamFamily(sCurrentLine.split("\\s+")[2]);
				} else if(sCurrentLine.startsWith("#=GS") && sCurrentLine.contains("AC")) { //new Prot
					currentFamily.addProtein(sCurrentLine.split("\\s+")[1].split("/")[0].intern());
				}
			}
			ret.add(currentFamily);
			br.close();
			
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(Global.PFAM_PATH+".bin"));
			oos.writeObject(ret);
			oos.close();
		}
		
		return ret;
	}

}
