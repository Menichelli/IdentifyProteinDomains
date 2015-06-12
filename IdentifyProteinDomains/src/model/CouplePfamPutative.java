/**
 * 
 */
package model;

/**
 * @author christophe
 *
 */
public class CouplePfamPutative {
	
	private PfamFamily pfam;
	private PutativeDomain putative;
	
	public CouplePfamPutative(PfamFamily pfam, PutativeDomain putative) {
		super();
		this.pfam = pfam;
		this.putative = putative;
	}
	
	public PfamFamily getPfam() {
		return pfam;
	}
	
	public void setPfam(PfamFamily pfam) {
		this.pfam = pfam;
	}
	
	public PutativeDomain getPutative() {
		return putative;
	}
	
	public void setPutative(PutativeDomain putative) {
		this.putative = putative;
	}
	
}
