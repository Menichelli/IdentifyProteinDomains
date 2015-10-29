package global;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import model.Proteome;
import tools.parser.RefProtParser;


public class Global {

	public static String PATHPROPERTIES;

	private static Global instance;

	//All
	public static boolean VERBOSE;
	public static boolean DYNAMIC_DISPLAY;
	public static boolean COMPUTE_FDR;
	public static String R_SCRIPT_PATH;
	public static Map<String, String> MAPPING_CLUSTERID;
	public static boolean KEEPONLYBESTPVALUE;
	public static Proteome PROTEOME_AIMED;
	public static boolean ONEMODELPERDOMAIN;
	public static boolean PRINTPUTATIVESFASTA;
	//BlastResultsParser
	public static String BLAST_RESULTS_PATH;
	//ProteomeParser
	public static String PROTEOME_PATH;
	//PfamParser
	public static String PFAM_PATH;
	//RefProtParser
	public static String REFPROT_PATH;
	//RParser
	public static String R_RESULTS_PFPT_PATH;
	public static String R_RESULTS_PTPT_PATH;
	//HitsGathering
	public static int SIZE_HITS_MIN;
	public static int NUMBER_OF_HITS_MIN;
	public static double OVERLAP_RATE_MIN;
	public static Set<String> EXCLUDE_SPECIES_SET;
	//StatsPrinter
	public static String STATS_PFPT_PATH;
	public static String STATS_PTPT_PATH;
	public static int NB_SEQ_INTERSECT;
	//FastaPrinter
	public static String FASTA_DIR;
	//ResultsPrinter
	public static String RESULTS_PATH;
	//FDREstimator
	public static int FDR_NB_REPEATS;
	public static String FDR_TMP_PATH;
	
	
	private Global() {}

	public static void init() {
		if(Global.instance==null) {
			Global.instance = new Global();
			Properties props = new Properties();
			try {
				props.load(Global.class.getClassLoader().getResourceAsStream(PATHPROPERTIES));

				Global.VERBOSE = Boolean.parseBoolean(props.getProperty("verbose"));
				Global.DYNAMIC_DISPLAY = Boolean.parseBoolean(props.getProperty("dynamic_display"));
				Global.COMPUTE_FDR = Boolean.parseBoolean(props.getProperty("compute_fdr"));
				Global.KEEPONLYBESTPVALUE = Boolean.parseBoolean(props.getProperty("keep_only_best_Pvalue"));
				Global.ONEMODELPERDOMAIN = Boolean.parseBoolean(props.getProperty("one_model_per_domain"));
				Global.PRINTPUTATIVESFASTA = Boolean.parseBoolean(props.getProperty("print_putatives_fasta"));
				Global.R_SCRIPT_PATH = props.getProperty("r_script_path");
				if(!new File(Global.R_SCRIPT_PATH).exists()) throw new IOException("R script file not found.");
				Global.BLAST_RESULTS_PATH = props.getProperty("blast_results_path");
				if(!new File(Global.BLAST_RESULTS_PATH).exists()) throw new IOException("BLAST results file not found.");
				Global.PROTEOME_PATH = props.getProperty("proteome_path");
				if(!new File(Global.PROTEOME_PATH).exists()) throw new IOException("Proteome file not found.");
				Global.PFAM_PATH = props.getProperty("pfam_path");
				if(!new File(Global.PFAM_PATH).exists()) throw new IOException("Pfam file not found.");
				Global.REFPROT_PATH = props.getProperty("refprot_path");
				if(!new File(Global.REFPROT_PATH).exists()) throw new IOException("RefProt file not found.");
				Global.SIZE_HITS_MIN = Integer.parseInt(props.getProperty("size_hits_min"));
				Global.NUMBER_OF_HITS_MIN = Integer.parseInt(props.getProperty("nb_hits_min"));
				Global.OVERLAP_RATE_MIN = Double.parseDouble(props.getProperty("overlap_rate_min"));
				Global.EXCLUDE_SPECIES_SET = new HashSet<String>();
				for(String s : props.getProperty("exclude_species").split("\\|")) {
					if(s!=null && !s.isEmpty())
					Global.EXCLUDE_SPECIES_SET.add(s);
				}
				Global.STATS_PFPT_PATH = props.getProperty("stats_pfam_path");
				Global.STATS_PTPT_PATH = props.getProperty("stats_nopfam_path");
				Global.NB_SEQ_INTERSECT = Integer.parseInt(props.getProperty("nb_seq_intersect"));
				Global.R_RESULTS_PFPT_PATH = props.getProperty("r_results_pfpt_path");
				Global.R_RESULTS_PTPT_PATH = props.getProperty("r_results_ptpt_path");
				Global.FASTA_DIR = props.getProperty("fasta_dir");
				Global.RESULTS_PATH = props.getProperty("results_path");
				Global.FDR_NB_REPEATS = Integer.parseInt(props.getProperty("fdr_nb_repeats"));
				Global.FDR_TMP_PATH = props.getProperty("fdr_tmp_path");
				Global.MAPPING_CLUSTERID = RefProtParser.getClusterProtMapping();
				
				if(Boolean.parseBoolean(props.getProperty("print_options"))) {
					System.out.println("***\nUsing options:\n");
					System.out.println("Min nb sequences in commom: "+NB_SEQ_INTERSECT);
					System.out.println("Min domain size : "+SIZE_HITS_MIN);
					System.out.println("Min nb hits : "+NUMBER_OF_HITS_MIN);
					System.out.println("Hit min cover rate with domain : "+OVERLAP_RATE_MIN);
					System.out.println("Keep only best Pvalue: "+KEEPONLYBESTPVALUE);
					System.out.println("One model per domain : "+ONEMODELPERDOMAIN);
					System.out.println("***");
				}
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
	}
}
