package model;

public class CouplePutativePutative {
	
	private PutativeDomain putative1;
	private PutativeDomain putative2;
	public CouplePutativePutative(PutativeDomain putative1,
			PutativeDomain putative2) {
		super();
		this.putative1 = putative1;
		this.putative2 = putative2;
	}
	
	public PutativeDomain getPutative1() {
		return putative1;
	}
	public void setPutative1(PutativeDomain putative1) {
		this.putative1 = putative1;
	}
	public PutativeDomain getPutative2() {
		return putative2;
	}
	public void setPutative2(PutativeDomain putative2) {
		this.putative2 = putative2;
	}
	
}
