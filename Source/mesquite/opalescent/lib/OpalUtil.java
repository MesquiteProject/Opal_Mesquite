/* Opalescent source code.  Copyright 2012 Travis Wheeler and David Maddison.
Version 2.10, January 2012.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
*/
package mesquite.opalescent.lib;

import mesquite.categ.lib.*;

public class OpalUtil {

	//ints for DNA data
	/*The following are the ints used inside opal for nucleotides/ambiguity-codes
	  They are based on the order in this string:
        "ACGTUKMRYSWBVHDN"
         0123456789111111
                   012345  */
	
	static final int opalDNAGap = -2;   
	static final int opalDNAMissing = 15;   
	static final int opalDNAA = 0;
	static final int opalDNAC = 1;
	static final int opalDNAG = 2;
	static final int opalDNAT = 3;
                       
	//IUPAC DNA ambiguity codes
	static final int opalDNAR = 7;
	static final int opalDNAY = 8;
	static final int opalDNAB = 11;
	static final int opalDNAD = 14;
	static final int opalDNAH = 13;
	static final int opalDNAK = 5;
	static final int opalDNAM = 6;
	static final int opalDNAS = 9;
	static final int opalDNAV = 12;
	static final int opalDNAW = 10;

	//ints for protein data
	/*The following are the ints used inside opal for amino acids
	  They are based on the order in this string:
      "ARNDCQEGHILKMFPSTWYVBZX"
       01234567891111111111222
                 0123456789012  */

	static final int opalProteinGap = -2;   
	static final int opalProteinMissing = 22;   
	static final int opalProteinALA = 0;  
	static final int opalProteinARG = 1;
	static final int opalProteinASN = 2;
	static final int opalProteinASP = 3; 
	static final int opalProteinCYS = 4; 
	static final int opalProteinGLN = 5; 
	static final int opalProteinGLU = 6; 
	static final int opalProteinGLY = 7; 
	static final int opalProteinHIS = 8; 
	static final int opalProteinILE = 9; 
	static final int opalProteinLEU = 10; 
	static final int opalProteinLYS = 11;
	static final int opalProteinMET = 12; 
	static final int opalProteinPHE = 13; 
	static final int opalProteinPRO = 14;
	static final int opalProteinSER = 15; 
	static final int opalProteinTHR = 16; 
	static final int opalProteinTRP = 17; 
	static final int opalProteinTYR = 18;
	static final int opalProteinVAL = 19; 
	static final int opalProteinASX = 20; //(ASN/ASP = N/D)
	static final int opalProteinGLX = 21; //(GLN/GLU = Q/E)
	
	/*.................................................................................................................*/
	/** converts the sequences as given in OPALs int[][] format into the long[][] Mesquite expects */
	public static long[][] getMesquiteArray(int[][] opalResults, boolean isProtein) {
		if (opalResults==null)
			return null;
		
		int[][] flipped = new int[opalResults[0].length][opalResults.length];
		long[][] array = new long[opalResults[0].length][opalResults.length];
		
		for (int i=0; i<opalResults.length; i++) 
			for (int j=0; j<opalResults[0].length; j++)
				flipped[j][i] = opalResults[i][j];
		
		for (int i=0; i<array.length; i++) {
			for (int base=0;base<array[i].length;base++)
				array[i][base] = CategoricalState.inapplicable;   //set them all to inapplicable to start
		}

		for (int i=0; i<flipped.length; i++) {
			for (int base=0;base<flipped[i].length;base++) {
				int opalStates = flipped[i][base];
				if (isProtein) {
					switch (opalStates) {
					case opalProteinGap : 
						array[i][base]= CategoricalState.inapplicable;
						break;
					case opalProteinMissing: 
						array[i][base]= CategoricalState.unassigned;
						break;
					case opalProteinALA : 
						array[i][base]= ProteinState.makeSet(ProteinData.ALA);   
						break;
					case opalProteinARG : 
						array[i][base]= ProteinState.makeSet(ProteinData.ARG);   
						break;
					case opalProteinASN : 
						array[i][base]= ProteinState.makeSet(ProteinData.ASN);   
						break;
					case opalProteinASP : 
						array[i][base]= ProteinState.makeSet(ProteinData.ASP);   
						break;
					case opalProteinCYS : 
						array[i][base]= ProteinState.makeSet(ProteinData.CYS);   
						break;
					case opalProteinGLN : 
						array[i][base]= ProteinState.makeSet(ProteinData.GLN);   
						break;
					case opalProteinGLU : 
						array[i][base]= ProteinState.makeSet(ProteinData.GLU);   
						break;
					case opalProteinGLY : 
						array[i][base]= ProteinState.makeSet(ProteinData.GLY);   
						break;
					case opalProteinHIS : 
						array[i][base]= ProteinState.makeSet(ProteinData.HIS);   
						break;
					case opalProteinILE : 
						array[i][base]= ProteinState.makeSet(ProteinData.ILEU);   
						break;
					case opalProteinLEU : 
						array[i][base]= ProteinState.makeSet(ProteinData.LEU);   
						break;
					case opalProteinLYS : 
						array[i][base]= ProteinState.makeSet(ProteinData.LYS);   
						break;
					case opalProteinMET : 
						array[i][base]= ProteinState.makeSet(ProteinData.MET);   
						break;
					case opalProteinPHE : 
						array[i][base]= ProteinState.makeSet(ProteinData.PHE);   
						break;
					case opalProteinPRO : 
						array[i][base]= ProteinState.makeSet(ProteinData.PRO);   
						break;
					case opalProteinSER : 
						array[i][base]= ProteinState.makeSet(ProteinData.SER);   
						break;
					case opalProteinTHR : 
						array[i][base]= ProteinState.makeSet(ProteinData.THR);   
						break;
					case opalProteinTRP : 
						array[i][base]= ProteinState.makeSet(ProteinData.TRP);   
						break;
					case opalProteinTYR : 
						array[i][base]= ProteinState.makeSet(ProteinData.TYR);   
						break;
					case opalProteinVAL : 
						array[i][base]= ProteinState.makeSet(ProteinData.VAL);   
						break;
					case opalProteinASX : 
						array[i][base]= ProteinState.makeSet(ProteinData.ASN, ProteinData.ASP);
						break;
					case opalProteinGLX : 
						array[i][base]= ProteinState.makeSet(ProteinData.GLN, ProteinData.GLU);   
						break;
					default:
						array[i][base]= CategoricalState.unassigned;

					}

				} else {  // it's DNA data
					switch (opalStates) {
					case opalDNAGap : 
						array[i][base]= CategoricalState.inapplicable;
						break;
					case opalDNAMissing: 
						array[i][base]= CategoricalState.unassigned;
						break;
					case opalDNAA : 
						array[i][base]= DNAState.A;
						break;
					case opalDNAC : 
						array[i][base]= DNAState.C;
						break;
					case opalDNAG : 
						array[i][base]= DNAState.G;
						break;
					case opalDNAT : 
						array[i][base]= DNAState.T;
						break;
					case opalDNAR : 
						array[i][base]= CategoricalState.setUncertainty(DNAState.A | DNAState.G, true);
						break;
					case opalDNAY:
						array[i][base]= CategoricalState.setUncertainty(DNAState.C | DNAState.T, true);
						break;
					case opalDNAB:
						array[i][base]= CategoricalState.setUncertainty(DNAState.C | DNAState.G | DNAState.T, true);
						break;
					case opalDNAD:
						array[i][base]= CategoricalState.setUncertainty(DNAState.A | DNAState.G | DNAState.T, true);
						break;
					case opalDNAH:
						array[i][base]= CategoricalState.setUncertainty(DNAState.A | DNAState.C | DNAState.T, true);
						break;
					case opalDNAK:
						array[i][base]= CategoricalState.setUncertainty(DNAState.G | DNAState.T, true);
						break;
					case opalDNAM:
						array[i][base]= CategoricalState.setUncertainty(DNAState.A | DNAState.C, true);
						break;
					case opalDNAS:
						array[i][base]= CategoricalState.setUncertainty(DNAState.C | DNAState.G, true);
						break;
					case opalDNAV:
						array[i][base]= CategoricalState.setUncertainty(DNAState.A | DNAState.C | DNAState.G , true);
						break;
					case opalDNAW:
						array[i][base]= CategoricalState.setUncertainty(DNAState.A | DNAState.T, true);
						break;

					default:
						array[i][base]= CategoricalState.unassigned;

					}
				}
			}
		}

		return array;
	}
	/*.................................................................................................................*/
	/** converts the sequences as given in Mesquite's long format into the int[][] OPAL expects */
	public static int[][] getOpalArray(long[][] sequences, boolean isProtein, int[] mesquiteRows) {
		/*since empty rows in the Meqsuite matrix will not be copied into the Opal matrix,
		 * mesquiteRows will be used to track the relationship between rows in each matrix  */
		
		if (sequences==null)
			return null;
//		int max = 0;
//		for (int i=0; i<sequences.length; i++)
//			if (max<sequences[i].length)
//				max=sequences[i].length;
		//int[][] array = new int[sequences.length][max];
		int[][] array =  new int[sequences[0].length][sequences.length];
		int state = -1;

		int goodRowCnt = 0;
		for (int seqID=0;seqID<sequences[0].length;seqID++) {
			boolean hasCharacter = false;
			for (int colId=0; colId<sequences.length; colId++) {
				long states = sequences[colId][seqID];
				if (isProtein) {
					if (CategoricalState.isUnassigned(states)) {
						hasCharacter = true;
						array[seqID][colId]=opalProteinMissing; 
					} else if (CategoricalState.isInapplicable(states)) {
						array[seqID][colId]=opalProteinGap; 
					} else if (CategoricalState.hasMultipleStates(states)) {   // it has multiple states
						hasCharacter = true;
						for (state=0;state<CategoricalState.getMaxPossibleStateStatic(); state++) {
							//Did I do this right??? 
							//  do I need to test that there are only two states set?  
							// Or can I rely on these two being the only possible settings?
							if (CategoricalState.isElement(states, ProteinData.ASN) && CategoricalState.isElement(states, ProteinData.ASP)) { 
								array[seqID][colId]= opalProteinASX;   
							} else if (CategoricalState.isElement(states, ProteinData.GLN) && CategoricalState.isElement(states, ProteinData.GLU)) {
								array[seqID][colId]= opalProteinGLX;
							}
						}
					} else {
						hasCharacter = true;
						state = CategoricalState.getOnlyElement(states);
						switch (state) {
						case ProteinData.ALA : 
							array[seqID][colId]= opalProteinALA;   
							break;
						case ProteinData.ARG : 
							array[seqID][colId]= opalProteinARG;   
							break;
						case ProteinData.ASN : 
							array[seqID][colId]= opalProteinASN;   
							break;
						case ProteinData.ASP : 
							array[seqID][colId]= opalProteinASP;   
							break;
						case ProteinData.CYS : 
							array[seqID][colId]= opalProteinCYS;   
							break;
						case ProteinData.GLN : 
							array[seqID][colId]= opalProteinGLN;   
							break;
						case ProteinData.GLU : 
							array[seqID][colId]= opalProteinGLU;   
							break;
						case ProteinData.GLY : 
							array[seqID][colId]= opalProteinGLY;   
							break;
						case ProteinData.HIS : 
							array[seqID][colId]= opalProteinHIS;   
							break;
						case ProteinData.ILEU : 
							array[seqID][colId]= opalProteinILE;   
							break;
						case ProteinData.LEU : 
							array[seqID][colId]= opalProteinLEU;   
							break;
						case ProteinData.LYS : 
							array[seqID][colId]= opalProteinLYS;   
							break;
						case ProteinData.MET : 
							array[seqID][colId]= opalProteinMET;   
							break;
						case ProteinData.PHE : 
							array[seqID][colId]= opalProteinPHE;   
							break;
						case ProteinData.PRO : 
							array[seqID][colId]= opalProteinPRO;   
							break;
						case ProteinData.SER : 
							array[seqID][colId]= opalProteinSER;   
							break;
						case ProteinData.THR : 
							array[seqID][colId]= opalProteinTHR;   
							break;
						case ProteinData.TRP : 
							array[seqID][colId]= opalProteinTRP;   
							break;
						case ProteinData.TYR : 
							array[seqID][colId]= opalProteinTYR;   
							break;
						case ProteinData.VAL : 
							array[seqID][colId]= opalProteinVAL;   
							break;
						default:
							array[seqID][colId]= opalProteinMissing;

						}
					}

				} else {  // it's DNA data
					if (states == DNAState.unassigned) {
						hasCharacter = true;
						array[seqID][colId]=opalDNAMissing; 
					} else if (states == DNAState.inapplicable) {
						array[seqID][colId]=opalDNAGap; 
					} else {
						hasCharacter = true;
						states = states & DNAState.statesBitsMask;
						if (states == DNAState.A)
							array[seqID][colId]=opalDNAA; 
						else if (states == DNAState.C)
							array[seqID][colId]=opalDNAC; 
						else if (states == DNAState.G)
							array[seqID][colId]=opalDNAG; 
						else if (states == DNAState.T)
							array[seqID][colId]=opalDNAT; 
						else if (states == (DNAState.A | DNAState.C))
							array[seqID][colId]=opalDNAM; 
						else if (states == (DNAState.A  | DNAState.G))  
							array[seqID][colId]=opalDNAR; 
						else if (states == (DNAState.A  | DNAState.T))  
							array[seqID][colId]=opalDNAW; 
						else if (states == (DNAState.C | DNAState.G ))  
							array[seqID][colId]=opalDNAS; 
						else if (states == (DNAState.C |  DNAState.T)  )
							array[seqID][colId]=opalDNAY; 
						else if (states == (DNAState.G | DNAState.T)  )
							array[seqID][colId]=opalDNAK; 
						else if (states == (DNAState.A | DNAState.C | DNAState.G )  )
							array[seqID][colId]=opalDNAV; 
						else if (states == (DNAState.A | DNAState.C | DNAState.T)  )
							array[seqID][colId]=opalDNAH; 
						else if (states == (DNAState.A |  DNAState.G | DNAState.T)  )
							array[seqID][colId]=opalDNAD; 
						else if (states == (DNAState.C | DNAState.G | DNAState.T) ) 
							array[seqID][colId]=opalDNAB; 
						else if (states == (DNAState.A | DNAState.C | DNAState.G | DNAState.T)  )
							array[seqID][colId]=opalDNAMissing;
					}
				}
			}
/*
			for (int base=sequences[0].length;base<max; base++) {
				if (isProtein)
					array[seqID][base] = opalProteinGap;
				else
					array[seqID][base] = opalDNAGap;
			}
	*/		
			if (hasCharacter) 
				mesquiteRows[goodRowCnt++] = seqID;
			
		}
		
		int[][] ret = new int[goodRowCnt][];
		for (int i=0; i<goodRowCnt; i++)
			ret[i] = array[mesquiteRows[i]];
		
		
		return ret;
	}

	public static int[][] addBackEmptyRows (int[][] opalResults, int[] mesquiteRows ) {
	
		int[][] returned = new int[mesquiteRows.length][];
		int row = 0;
		for (int i=0; i<mesquiteRows.length; i++){
			if (i==mesquiteRows[row]) {
				returned[i] = opalResults[row++];
			} else {
				//stick a row in with all gap characters
				returned[i] = new int[opalResults[0].length];
				for (int j=0; j<opalResults[0].length; j++)
					returned[i][j]=opalDNAMissing; 
			}
		}
		return returned;
	}



}
