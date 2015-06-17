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
public class ConservationStatsPrinter {

	private static ConservationStatsPrinter instance;
	private BufferedWriter writer;
	private String path;
	
	private ConservationStatsPrinter(String path) throws IOException {
		this.path = path;
		File f = new File(path);
		if(!f.exists()) f.createNewFile();
		FileWriter fw = new FileWriter(f);
		writer = new BufferedWriter(fw);
		writer.write("nbRetained nbTotal ratio\n");
	}

	public static ConservationStatsPrinter getInstance(String path) throws Exception {
		if(instance!=null && !path.equals(instance.path)) throw new Exception("You didn't close previous ConservationStatsPrinter!");
		instance = (instance==null)? new ConservationStatsPrinter(path):instance;
		return instance;
	}

	public void close() throws Exception {
		writer.close();
		instance = null;
	}

	public synchronized void addEntry(int nbRetained, int nbTotal) throws Exception {
		writer.append(nbRetained+" "+nbTotal+" "+(double)((double)nbRetained/(double)nbTotal)+"\n");
	}

}
