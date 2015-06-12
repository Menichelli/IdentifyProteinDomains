/**
 * 
 */
package model;

/**
 * @author christophe
 *
 */
public class BlastHit {
	private String queryName;
	private String querySpecies;
	private String subjectName;
	private String subjectSpecies;
	private double identity;
	private int alignmentLength;
	private int qStart;
	private int qEnd;
	private int sStart;
	private int sEnd;
	private double eValue;
	
	public BlastHit(String queryName, String querySpecies, String subjectName,
			String subjectSpecies, double identityPercent,
			int alignmentLength, int qStart, int qEnd,
			int sStart, int sEnd, double eValue) {
		super();
		this.queryName = queryName;
		this.querySpecies = querySpecies;
		this.subjectName = subjectName;
		this.subjectSpecies = subjectSpecies;
		this.identity = identityPercent;
		this.alignmentLength = alignmentLength;
		this.qStart = qStart;
		this.qEnd = qEnd;
		this.sStart = sStart;
		this.sEnd = sEnd;
		this.eValue = eValue;
	}
	/**
	 * @return the queryName
	 */
	public String getQueryName() {
		return queryName;
	}
	/**
	 * @param queryName the queryName to set
	 */
	public void setQueryName(String queryName) {
		this.queryName = queryName;
	}
	/**
	 * @return the querySpecies
	 */
	public String getQuerySpecies() {
		return querySpecies;
	}
	/**
	 * @param querySpecies the querySpecies to set
	 */
	public void setQuerySpecies(String querySpecies) {
		this.querySpecies = querySpecies;
	}
	/**
	 * @return the subjectName
	 */
	public String getSubjectName() {
		return subjectName;
	}
	/**
	 * @param subjectName the subjectName to set
	 */
	public void setSubjectName(String subjectName) {
		this.subjectName = subjectName;
	}
	/**
	 * @return the subjectSpecies
	 */
	public String getSubjectSpecies() {
		return subjectSpecies;
	}
	/**
	 * @param subjectSpecies the subjectSpecies to set
	 */
	public void setSubjectSpecies(String subjectSpecies) {
		this.subjectSpecies = subjectSpecies;
	}
	/**
	 * @return the identityPercent
	 */
	public double getIdentityPercent() {
		return identity;
	}
	/**
	 * @param identityPercent the identityPercent to set
	 */
	public void setIdentityPercent(double identityPercent) {
		this.identity = identityPercent;
	}
	/**
	 * @return the alignmentLength
	 */
	public int getAlignmentLength() {
		return alignmentLength;
	}
	/**
	 * @param alignmentLength the alignmentLength to set
	 */
	public void setAlignmentLength(int alignmentLength) {
		this.alignmentLength = alignmentLength;
	}
	/**
	 * @return the qStart
	 */
	public int getqStart() {
		return qStart;
	}
	/**
	 * @param qStart the qStart to set
	 */
	public void setqStart(int qStart) {
		this.qStart = qStart;
	}
	/**
	 * @return the qEnd
	 */
	public int getqEnd() {
		return qEnd;
	}
	/**
	 * @param qEnd the qEnd to set
	 */
	public void setqEnd(int qEnd) {
		this.qEnd = qEnd;
	}
	/**
	 * @return the sStart
	 */
	public int getsStart() {
		return sStart;
	}
	/**
	 * @param sStart the sStart to set
	 */
	public void setsStart(int sStart) {
		this.sStart = sStart;
	}
	/**
	 * @return the sEnd
	 */
	public int getsEnd() {
		return sEnd;
	}
	/**
	 * @param sEnd the sEnd to set
	 */
	public void setsEnd(int sEnd) {
		this.sEnd = sEnd;
	}
	/**
	 * @return the eValue
	 */
	public double geteValue() {
		return eValue;
	}
	/**
	 * @param eValue the eValue to set
	 */
	public void seteValue(double eValue) {
		this.eValue = eValue;
	}
	
	@Override
	public String toString() {
		return String
				.format("BlastHit [queryName=%s, querySpecies=%s, subjectName=%s, subjectSpecies=%s, identityPercent=%s, alignmentLength=%s, qStart=%s, qEnd=%s, sStart=%s, sEnd=%s, eValue=%s]",
						queryName, querySpecies, subjectName, subjectSpecies,
						identity, alignmentLength, qStart,
						qEnd, sStart, sEnd, eValue);
	}
	
}
