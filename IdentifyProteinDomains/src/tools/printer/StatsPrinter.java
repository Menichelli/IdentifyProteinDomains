/**
 * 
 */
package tools.printer;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * @author christophe
 *
 */
public class StatsPrinter {

	private static StatsPrinter instance;
	private BufferedWriter writer;
	private String path;
	
	private StatsPrinter(String path) throws IOException {
		this.path = path;
		File f = new File(path);
		if(!f.exists()) f.createNewFile();
		FileWriter fw = new FileWriter(f);
		writer = new BufferedWriter(fw);
		writer.write("domainName1 domainName2 nbProtDomain1 nbProtDomain2 nbProtIntersec\n");
	}

	public static StatsPrinter getInstance(String path) throws Exception {
		if(instance!=null && !path.equals(instance.path)) throw new Exception("You didn't close previous StatsPrinter!");
		instance = (instance==null)? new StatsPrinter(path):instance;
		return instance;
	}

	public void close() throws Exception {
		writer.close();
		instance = null;
	}

	public synchronized void addEntry(String domainName1, String domainName2, int nbProtName1, int nbProtDomain2, int nbProtIntersec) throws Exception {
		writer.append(domainName1+" "+domainName2+" "+nbProtName1+" "+nbProtDomain2+" "+nbProtIntersec+"\n");
	}

}
