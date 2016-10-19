/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.helpers;

import de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.samhelpers.SamSplitter;


/**
 * @author Jhoser
 *
 */
public class Sequence
{
	public static final int UNKNOWN 	= 0;
	public static final int WITHQUAL	= 1;
	public static final int DNA 		= 2;
	public static final int RNA			= 4;
	public static final int AA			= 8;
	private static final String NL = System.getProperty("line.separator");
	/**
	 * The id of the sequence (assesionnumber, chromosome...)
	 */
	private String id;
	/**
	 * the sequencestring itself
	 */
	private char[] sequence;
	
	/**
	 * Type of the Sequence with type being one of the 
	 * defined in Sequence.DNA Sequence.RNA Sequence.AA
	 * or Sequence.UNKNOWN (or any bitwise combination of those)
	 */
	private int seqType;
	
	/**
	 * if given contains the quality values for the sequence
	 */
	private char[] quality;

	/**
	 * Record if we are Zero (0) or One(1) based. For sequence offset.
	 */
	private boolean isZeroBased=true;
	
	/**
	 * Create a new Sequence-Object, using sequence ID only
	 * (allows later addition of sequence itself (and qual)
	 * @param id (fasta-)id of the sequence
	 */
	public Sequence(String id,boolean trimToFirst)
		{
			if(!trimToFirst)
				this.id = id;
			else
				{
					if(id.charAt(0)=='>')
						{
							id = id.substring(1);
						}
					int posOfSpace = id.indexOf(' ');
					if(posOfSpace>0)
						id=id.substring(0,posOfSpace);
					this.id = id;
				}
		}
	
	/**
	 * Create a new Sequence-Object, using sequence ID and sequence
	 * (allows later addition of qual - if given.)
	 * @param id (fasta-)id of the sequence
	 * @param sequence the sequence of this sequence
	 */
	public Sequence(String id, String sequence)
		{
			this.id = id;
			this.sequence = sequence.toCharArray();
			sequence = null;
		}
	
	/**
	 * Create a new Sequence-Object, using sequence ID, sequence and qualities (ascii encoded)
	 * 
	 * @param id (fasta-)id of the sequence
	 * @param sequence the sequence of this sequence
	 * @param qual the quality string of this sequence (ascii encoded quals)
	 * @throws ParameterInvalidException 
	 */
	public Sequence(String id, String sequence, String qual) throws ParameterInvalidException
		{
			this.id = id;
			this.sequence = sequence.toCharArray();
			this.quality = qual.toCharArray();
			if(this.sequence.length!=this.quality.length)
				{
					System.out.flush();
					System.err.flush();
					System.err.println("Given Sequence length and given quality length differ... cannot continue. Aborting.");
					if(!SamSplitter.noExit) 
						System.exit(2);
					else
						throw new IllegalArgumentException("Given Sequence length and given quality length differ... cannot continue. Aborting.");
				}
			sequence = null;
			qual = null;
		}
	
	/**
	 * Set the sequence itself
	 * @param seq the sequenceString to be stored
	 */
	public void setSequence(String seq)
		{
			this.sequence = seq.toCharArray();
		}
	
	/**
	 * Set the qualitySequence itself
	 * @param qual the qualitySequenceString to be stored
	 */
	public void setQualSequence(String qual)
		{
			this.quality = qual.toCharArray();
		}
	
	/**
	 * Set the sequenceType as defined via Sequence.STATICS
	 * (bitwise combinations are ok, to: e.g. 3 for DNA+Qual, 8 for AA only, 9 for AA+qual etc)
	 * @param type the type of the herin stored sequence
	 */
	public void setSequenceType(int type)
		{
			this.seqType = type;
		}
	
	/**
	 * Get the stored sequenceType, if none stored returns Sequence.UNKNOWN
	 * @return the sequenceType
	 */
	public int getSequenceType()
		{
			return this.seqType;
		}
	
	/**
	 * Returns the entire sequence in one String object 
	 * (for outputting etc)
	 * @return the entire Sequence
	 */
	public String getWholeSequence()
		{
			return new String(this.sequence);
		}
	
	/**
	 * Returns the entire sequence in one char-array object 
	 * (for outputting etc)
	 * @return the entire Sequence
	 */
	public char[] getWholeSequenceCharArr()
		{
			return this.sequence;
		}
	
	/**
	 * Returns all qualities as String (if there are any, otherwise returns null)
	 * @return String of qualities or null
	 */
	public String getWholeQuality()
		{
			if(this.quality.length>0)
				return new String(this.quality);
			else
				return null;
		}
	
	/**
	 * Returns all qualities as char-Arrray (if there are any, otherwise returns null)
	 * @return char-array of qualities or null
	 */
	public char[] getWholeQualityCharArr()
		{
			if(this.quality.length>0)
				return this.quality;
			else
				return null;
		}
	
	/**
	 * Get a single char from the sequence, the one whose position has been specified via pos
	 * @param pos the position from where in the sequence data (char) shall be fetched 
	 * @return the char found at sequenceposition[pos]
	 */
	public char getSequencePos(int pos)
		{
			if(isZeroBased)
				return this.sequence[pos];
			else
				{
					return this.sequence[pos-1];
				}
		}
	/**
	 * Get a single char from the qualitysequence, the one whose position has been specified via pos
	 * @param pos the position from where in the qualitysequence data (char) shall be fetched 
	 * @return the char found at qualitysequenceposition[pos]
	 */
	public char getQualityPos(int pos)
		{
			return this.quality[pos];
		}
	
	/**
	 * Method returns the overall length of the herin stored sequence.
	 * @return the length of the sequence (and if set of the qualities)
	 */
	public int getLen()
		{
			return this.sequence.length;
		}
	
	/**
	 * Returns the id of this sequence
	 * @return the id of this sequence
	 */
	public String getId()
		{
			return this.id;
		}
	
	/**
	 * Method allows requesting a part of the sequence.
	 * from is the inclusive start (0-based count)
	 * to is the inclusive end (0-based count)
	 * @param from startpoint (included) from where to begin taking the sequence (0-based)
	 * @param to endpoint (included) to what position we want the sequence. (0-based)
	 * @return String containing the requested subsequence or null if the request was to large/invalid
	 */
	public String getSequencePart(int from, int to)
		{
			if(from>=0 && to <= this.sequence.length)
				{
					char[] out = new char[to-from+1];
					int p=0;
					if(isZeroBased)
						{
							for(int i=from;i<=to;i++)
								{
									out[p] = this.sequence[i];
									//p +=1;
									++p;
								}
						}
					else
						{
							for(int i=from-1;i<=to-1;i++)
								{
									out[p] = this.sequence[i];
									//p +=1;
									++p;
								}
						}
					return new String(out);
				}
			else
				{
					System.err.println("Request for subsequence from pos.:" + from + " to pos.:" + to + " is invalid for sequence of length " + this.sequence.length);
					return null;
				}
		}
	
	/**
	 * Method allows requesting the complementary part of the sequence.
	 * from is the inclusive start (0-based count) in the orig. sequence
	 * to is the inclusive end (0-based count) in the orig. sequence
	 * @param from startpoint (included) from where to begin taking the comp.-sequence (0-based)
	 * @param to endpoint (included) to what position we want the comp.-sequence. (0-based)
	 * @return String containing the requested subsequence or null if the request was to large/invalid
	 */
	public String getComplementarySequencePart(int from, int to)
		{
			if(from>=0 && to <= this.sequence.length)
				{
					char[] out = new char[to-from+1];
					int p=0;
					for(int i=to;i>=from;i--)
						{
							out[p] = this.getInverseNT(this.sequence[i]);
							p +=1;
						}
					return new String(out);
				}
			else
				{
					System.err.println("Request for subsequence from pos.:" + from + " to pos.:" + to + " is invalid for sequence of length " + this.sequence.length);
					return null;
				}
		}
	
	/**
	 * Compute the complementary Sequence to our already stored sequence
	 * and return it. (as String)
	 * @return complementary Sequence (DNA) to stored sequence.
	 */
	public String getWholeComplementarySequence()
		{
			int len = this.sequence.length;
			char[] compSeq = new char[len];
			for(int i=0;i<len;i++)
				{
					compSeq[i] = getInverseNT(this.sequence[len-1-i]);
				}
			return new String(compSeq);
			
		}
	
	/**
     * 'Compute' inverse NT
     * @param nt the Nucleotide to invert..
     * @return its nt-binding partner
     */
    private char getInverseNT(char nt)
        {
            if(nt=='C')
                return 'G';
            else if(nt=='G')
                return 'C';
            else if(nt=='A')
                return 'T';
            else if(nt=='T')
                return 'A';
            else
                return nt;
        }

	/**
	 * @return the isZeroBased
	 */
	public boolean isZeroBased()
		{
			return isZeroBased;
		}

	/**
	 * Setting this to true means that a position-call for [0] will return the first base of the sequence,
	 * setting it to false, means calling [0] will fail, and [1] will return the first base of the sequence.
	 * @param isZeroBased the isZeroBased to set
	 */
	public void setIsZeroBased(boolean isZeroBased)
		{
			this.isZeroBased = isZeroBased;
		}

	/**
	 * Returns this Sequence' Objects properties as a Fasta/Fastq formatted string
	 * Depends on quality present, or not, and of course on the flag:
	 * @param printQual Whether or not to print output with quality or without: backDoor Fastq/fasta converter
	 */
	public String toString(boolean printQual)
		{
			boolean pQ = false;
			if(printQual)
				if(this.quality==null)
					{
						System.err.println("Cannot print quality: No quality values are stored!");
					}
				else
					{
						pQ = true;
					}
			
			StringBuilder sb = new StringBuilder();
			if(pQ)
				{
					if(this.id.charAt(0)!='@')
						sb.append('@');
				}
			else
				{
					if(this.id.charAt(0)!='>')
						sb.append('>');
				}
			
			sb.append(this.id);
			sb.append(Sequence.NL);
			sb.append(this.sequence);
			if(pQ)
				{
					sb.append(Sequence.NL);
					sb.append('+');
					sb.append(Sequence.NL);
					sb.append(this.quality);
				}
			
			return sb.toString();
		}
	
	
	
	/**
	 * Method to trim this sequence' Objects sequence to 1 nt!
	 * For keeping PE-information!
	 */
	public void trimSequenceRegionTo1nt()
		{
			if(this.sequence.length>0)
				{
					if(this.isZeroBased)
						this.trimSequenceRegion(1, this.sequence.length);
					else
						this.trimSequenceRegion(2, this.sequence.length);
				}
			else
				{
					this.sequence = new char[1];
					this.sequence[0] = 'N';
					this.quality = new char[1];
					this.quality[0] = '!';
				}
		}
	
	/**
	 * Method to trim this sequence' Objects sequence.
	 * Parameters specifiy which region to trim/erase/delete
	 * trimSequenceRegion(0,4) does e.g. remove the first 4 nt/aa/... if we are zero-based,
	 * @param fromPos
	 * @param toPos
	 */
	public void trimSequenceRegion(int fromPos, int toPos)
		{
			if(fromPos>=0 && toPos <= this.sequence.length)
				{
					//Seq is 20bp
					//trim 0,4 -> copy 4...this.sequence.length=20
					//char[] newseq = Arrays.copyOfRange(this.sequence, from, to)
					int newLen = this.sequence.length - (toPos-fromPos);
					boolean wQual = (this.quality!=null);
					
					char[] newSeq = new char[newLen];
					char[] newQual = null;
					if(wQual)
						newQual = new char[newLen];
					boolean isAtStart;
					
					if(this.isZeroBased)
						isAtStart = (fromPos==0);
					else
						isAtStart = (fromPos==1);
					
					if(isAtStart)
						{
							System.arraycopy(this.sequence, toPos, newSeq, 0, newLen);
							if(wQual)
								System.arraycopy(this.quality, toPos, newQual, 0, newLen);
						}
					else
						{
							System.arraycopy(this.sequence, 0, newSeq, 0, newLen);
							if(wQual)
								System.arraycopy(this.quality, 0, newQual, 0, newLen);
						}
					this.sequence = newSeq;
					this.quality = newQual;
				}
			else 
				{
					if(this.sequence.length==0)
						return;
					
					System.out.println("Sequence.trim: seq: " + (this.isZeroBased?"0":"1") + " - " + this.sequence.length + "(len) \t request: " + fromPos + " - " + toPos);
					//System.err.println("Not trimming, from and to are not valid (ie. beyond sequence range)");
					
					if(fromPos>this.sequence.length)
						{
							return;
						}
					if(fromPos<0)
						return;
					
					
					if(toPos>=this.sequence.length)
						{
							this.trimSequenceRegion(fromPos, this.sequence.length-1);
						}
				}
		}
	
	/**
	 * Trim the sequence starting from the position where a char (as given to the function)
	 * is encountered. Trim that and any position beyond that from the read.
	 * @param bchar the char to look for.
	 * @return the number of positions that were trimmed.
	 */
	public int trimQualB(char bchar)
		{
			for(int i=0;i<this.quality.length;++i)
				{
					if(this.quality[i]!=bchar)
						{
							continue;
						}
					else
						{
							//found B. discard anything beyond.
							this.trimSequenceRegion(i, this.quality.length);
							return (this.quality.length - i);
						}
				}
			return 0;
		}

	/**
	 * Trim the sequence starting from the position where a quality (as given to the function)
	 * is encountered. Since we have an encoding in place, we need to decode it using the given offset.
	 * Trim that and any position beyond that from the read.
	 * @param bqual the quality to look for.
	 * @param offset the offset for the decoding.
	 * @return the number of positions that were trimmed.
	 */
	public int trimQualValue(int bqual, int offset)
		{
			for(int i=0;i<this.quality.length;++i)
				{
					if((this.quality[i]-offset)!=bqual)
						{
							continue;
						}
					else
						{
							//found B. discard anything beyond.
							this.trimSequenceRegion(i, this.quality.length);
							return (this.quality.length - i);
						}
				}
			return 0;
		}
	
	public int trimPolyAT(boolean trimStart, boolean trimEnd, int minPolyATsize)
		{
			int startPos = 0;
			int endPos = this.sequence.length;
			int trimCounter =0;
			char trimChar = 1;
			
			if(trimStart)
				{
					if(this.sequence[0]=='A' || this.sequence[0]=='T')
						{
							trimChar = this.sequence[0];
							++trimCounter;
							while(true)
								{
									if(this.sequence[trimCounter]==trimChar)
										{
											++trimCounter;
										}
									else
										{
											if(trimCounter>=minPolyATsize)
												startPos = trimCounter;
											trimCounter = 0;
											trimChar=1;
											break;
										}
								}
						}
				}
			
			if(trimEnd)
				{
					int i=1;
					if(this.sequence[endPos-i]=='A' || this.sequence[endPos-i]=='T')
						{
							trimChar = this.sequence[endPos-i];
							++trimCounter;
							while(true)
								{
									if((endPos-(trimCounter+i))>=0 && this.sequence[endPos-(trimCounter+i)]==trimChar)
										{
											++trimCounter;
										}
									else
										{
											//TODO:skip lowqual - other than AT char. and continue.
											if(trimCounter>=minPolyATsize)
												endPos = endPos-trimCounter;
											trimCounter = 0;
											trimChar=1;
											break;
										}
								}
						}
				}
			
			int sumTrim = startPos+(this.sequence.length-endPos);
			
			//System.out.println("Trimming:");
			//System.out.println(this.toString(false));
			//System.out.println("to");
			
			if(startPos>0)
				this.trimSequenceRegion(0, startPos);
			if(endPos<this.sequence.length)
				this.trimSequenceRegion(endPos, this.sequence.length);
			
			//this.trimSequenceRegion(startPos, endPos);
			//System.out.println(this.toString(false));
			
			return sumTrim;
		}
	
	public static void main(String[] args) throws ParameterInvalidException
		{
			String id="@test1";
			String seq1="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
			String qul1="hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhMMMMMMMMMMMMMMMM";
			Sequence s = new Sequence(id,seq1,qul1);
			System.out.println(s.toString(true));
			
			Sequence s2 = new Sequence(id,seq1,qul1);
			s2.setIsZeroBased(true);
			System.out.println(s2.toString(true));
			
			s.trimSequenceRegion(0, 5);
			System.out.println("Length: " + s.getLen());
			System.out.println(s.toString(true));
			
			s2.trimSequenceRegion(40, 48);
			System.out.println("Length: " + s2.getLen());
			System.out.println(s2.toString(true));
			
			System.out.println("Test:");
			s2.trimSequenceRegionTo1nt();
			System.out.println("Length: " + s2.getLen());
			System.out.println(s2.toString(true));
			
			s.trimSequenceRegionTo1nt();
			System.out.println("Length: " + s.getLen());
			System.out.println(s.toString(true));
		}

	
}

