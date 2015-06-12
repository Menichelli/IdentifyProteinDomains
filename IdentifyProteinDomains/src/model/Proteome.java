package model;

import java.util.HashSet;
import java.util.Set;

public class Proteome {
	
	Set<Protein> proteins;
	
	public Proteome() {
		super();
		this.proteins = new HashSet<Protein>();
		
	}
	
	public Set<Protein> getProteins() {
		return this.proteins;
	}
	
	public void addProtein(Protein p) {
		this.proteins.add(p);
	}
	
	public void removeProtein(Protein p) {
		this.proteins.remove(p);
	}
	
	public Protein getProteinByName(String protName) {
		Protein ret = null;
		for(Protein p : this.proteins) {
			if((p.getProtName()+"_"+p.getSpecies()).equals(protName)) {ret = p; break;}
		}
		return ret;
	}
	
}
