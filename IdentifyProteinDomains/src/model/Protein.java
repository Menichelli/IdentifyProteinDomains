package model;

public class Protein {
	
	private String protName;
	private String protIdentifier;
	private String species;
	private String sequence;
	
	/**
	 * @param protName
	 * @param protIdentifier
	 * @param species
	 * @param sequence
	 */
	public Protein(String protName, String protIdentifier, String species,
			String sequence) {
		super();
		this.protName = protName;
		this.protIdentifier = protIdentifier;
		this.species = species;
		this.sequence = sequence;
	}

	public String getProtName() {
		return protName;
	}

	public void setProtName(String protName) {
		this.protName = protName;
	}

	public String getProtIdentifier() {
		return protIdentifier;
	}

	public void setProtIdentifier(String protIdentifier) {
		this.protIdentifier = protIdentifier;
	}

	public String getSpecies() {
		return species;
	}

	public void setSpecies(String species) {
		this.species = species;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	public String toString() {
		return this.protName;
	}
	
}
