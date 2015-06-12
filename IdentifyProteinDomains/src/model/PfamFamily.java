package model;

import java.io.Serializable;
import java.util.HashSet;
import java.util.Set;

public class PfamFamily implements Serializable {

	private static final long serialVersionUID = 1L;
	private String familyName;
	private Set<String> proteins; //Store only the set of prot names

	public PfamFamily(String familyName) {
		this.familyName = familyName;
		this.proteins = new HashSet<String>();
	}

	public void addProtein(String protName) {
		this.proteins.add(protName);
	}

	public Set<String> getAllProteinNames() {
		return this.proteins;
	}

	public String getFamilyName() {
		return this.familyName;
	}

	public boolean isFoundOn(String protName) {
		boolean ret=false;
		ret = proteins.contains(protName);
		return ret;
	}

}
