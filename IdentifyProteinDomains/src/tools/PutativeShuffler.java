package tools;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import model.BlastHit;
import model.PutativeDomain;

public class PutativeShuffler {

	private static PutativeShuffler instance;
	private Random random;

	private PutativeShuffler() {
		random = new Random();
	}

	public static PutativeShuffler getInstance() {
		instance = (instance==null)?new PutativeShuffler():instance;
		return instance;
	}

	public Map<String,Set<PutativeDomain>> shuffle(final Map<String,Set<PutativeDomain>> putativeDomainsByProt) {
		Map<String,Set<PutativeDomain>> ret = new HashMap<String, Set<PutativeDomain>>(putativeDomainsByProt);

		String[] protNames = (String[])ret.keySet().toArray(new String[0]);
		Set<BlastHit> tmp;
		int rdmProt,rdmD;
		Set<PutativeDomain> listDomains;
		PutativeDomain d;
		for(int index = 0; index < protNames.length-1; index++) {
			for(PutativeDomain domain : ret.get(protNames[index])) {
				rdmProt = pickRandom(index+1, protNames.length-1);
				listDomains = ret.get(protNames[rdmProt]);
				rdmD = pickRandom(0, listDomains.size()-1);
				d = ((PutativeDomain[])listDomains.toArray(new PutativeDomain[0]))[rdmD];
				tmp = domain.getBlastHits();
				domain.setBlastHits(d.getBlastHits());
				d.setBlastHits(tmp);
			}
		}
		return ret;
	}

	private int pickRandom(int min, int max) {
		return random.nextInt(max+1-min)+min;
	}

}
