package model;

public class ValidatedDomain {
	
	private String identifierValidatedDomain;
	private String identifierValidatingDomain;
	private double score;
	
	public ValidatedDomain(String identifierValidatedDomain,
			String identifierValidatingDomain, double score) {
		super();
		this.identifierValidatedDomain = identifierValidatedDomain;
		this.identifierValidatingDomain = identifierValidatingDomain;
		this.score = score;
	}
	
	public String getIdentifierValidatedDomain() {
		return identifierValidatedDomain;
	}
	public void setIdentifierValidatedDomain(String identifierValidatedDomain) {
		this.identifierValidatedDomain = identifierValidatedDomain;
	}
	public String getIdentifierValidatingDomain() {
		return identifierValidatingDomain;
	}
	public void setIdentifierValidatingDomain(String identifierValidatingDomain) {
		this.identifierValidatingDomain = identifierValidatingDomain;
	}
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
	}

}
