/**
 * 
 */
package de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.assembly;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import de.helmholtz_muenchen.ibis.ThreadedReadFilter.utils.FileHelpers;
import de.helmholtz_muenchen.ibis.ThreadedReadFilter.utils.ParameterInvalidException;
import de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.helpers.Sequence;


/**
 * @author jonathan.hoser
 *
 */
public class TIGRrepairer {
	private static final String OUTDIR_STR="--outDir=";
	private static final String INFILE_STR="--in=";
	private File inFile;
	private File outFileForward;
	private File outFileReverse;
	private File outFileSingle;

	
	private void setInFile(String strInFile) throws ParameterInvalidException 
		{
			File _inFile = FileHelpers.checkInFile(strInFile);
			this.inFile = _inFile;
		}

	private void setOutDir(String strOutDir) throws ParameterInvalidException 
		{
			if(strOutDir!=null)
				{
					File _outDir = FileHelpers.checkOutDir(strOutDir);
					
					this.outFileReverse = new File(_outDir.getAbsolutePath()+System.getProperty("file.separator")+this.inFile.getName()+".rev");
					this.outFileForward = new File(_outDir.getAbsolutePath()+System.getProperty("file.separator")+this.inFile.getName()+".fwd");
					this.outFileSingle = new File(_outDir.getAbsolutePath()+System.getProperty("file.separator")+this.inFile.getName()+".sngl");
				}
			else
				{
					if(this.inFile == null)
						throw new ParameterInvalidException("Error 1: inFile CANNOT be empty! But it is");
					
					//grab the parent-File/Dir for OutDir.
					
					this.outFileReverse = new File(this.inFile.getParentFile().getAbsolutePath()+System.getProperty("file.separator")+this.inFile.getName()+".rev");
					this.outFileForward = new File(this.inFile.getParentFile().getAbsolutePath()+System.getProperty("file.separator")+this.inFile.getName()+".fwd");
					this.outFileSingle = new File(this.inFile.getParentFile().getAbsolutePath()+System.getProperty("file.separator")+this.inFile.getName()+".sngl");
				}
			System.out.println("Writing to:");
			System.out.println("FWD: " + this.outFileForward.getAbsolutePath());
			System.out.println("REV: " + this.outFileReverse.getAbsolutePath());
		}
	
	/**
	 * Workhorse class. Reads the file,
	 * gets the sequences, pairs them, and outputs them in an orderly fashion.
	 */
	private void readAndRePair() 
		{
			String line;
			long a = System.nanoTime();
			
			BufferedWriter bwFwd,bwRev,bwSngl;
			BufferedReader br;
			
			int curPos=0,newPos=0;
			
			Sequence s;
			String tigrStr,deciderStr;
			
			try
				{
		//			@2007-12-17_201B2AAXX_celegans_2_1_435_274
		//			GAAATTCACTACTATTAACAGAAGGACACTCATTAG
		//			+
		//			>>>>>>>>>>>>>>><;;;;;8887773.2&*-'+'
					
					br = new BufferedReader(new FileReader(this.inFile));
					bwFwd = new BufferedWriter(new FileWriter(this.outFileForward));
					bwRev = new BufferedWriter(new FileWriter(this.outFileReverse));
					
					HashMap<String,Sequence> fwd = new HashMap<String,Sequence>();
					HashMap<String,Sequence> rev = new HashMap<String,Sequence>();
					
					//char q = '\n';
					while((line=br.readLine())!=null) //name line
						{
							if(line.charAt(0)=='@')
								{
									//ID-line:
									s = new Sequence(line, false);
									//read sequence-line
									line=br.readLine();
									s.setSequence(line);
									//read extra line.
									line=br.readLine();
									//read qual-line.
									line=br.readLine();
									s.setQualSequence(line);
									
									//extract TIGR-id:
									//@TIGR_C1BA001TR_926918_1103683007942
									if(s.getId().startsWith("@TIGR"))
										{
											curPos = s.getId().indexOf('_');
											newPos = s.getId().indexOf('_', curPos+1);
											tigrStr = s.getId().substring(curPos, newPos-2);
											deciderStr = s.getId().substring(newPos-2, newPos);
											if(deciderStr.charAt(1)=='F')
												{
													if(rev.containsKey(tigrStr))
														{
															//write out!
															bwFwd.write(s.toString(true));
															bwFwd.newLine();
															s=rev.remove(tigrStr);
															bwRev.write(s.toString(true));
															bwRev.newLine();
														}
													else
														fwd.put(tigrStr, s);
												}
											if(deciderStr.charAt(1)=='R')
												{
													if(fwd.containsKey(tigrStr))
														{
															//write out!
															bwRev.write(s.toString(true));
															bwRev.newLine();
															s=fwd.remove(tigrStr);
															bwFwd.write(s.toString(true));
															bwFwd.newLine();
														}
													else
														rev.put(tigrStr, s);
												}
										}
									else
										{
											System.err.println("Sorry, we only do TIGR for now.");
											System.err.println("Encountered sequence ID = >"+s.getId()+"<");
											br.close();
											bwFwd.close();
											bwRev.close();
											System.exit(3);
										}
								}
						}
					
					bwSngl = new BufferedWriter(new FileWriter(this.outFileSingle));
					int snglCounter =0;
					if(!fwd.isEmpty())
						{
							System.err.println("FWD is not empty! Something is BAD, BAD, BAD!");
							ArrayList<String> list = new ArrayList<String>();
							list.addAll(fwd.keySet());
							Collections.sort(list);
							for (String id:list)
								{
									s = fwd.get(id);
									bwSngl.write(s.toString(true));
									bwSngl.newLine();
									++snglCounter;
									System.out.println(s.getId());
								}
							System.out.println("#####\nEntries: "+fwd.size()+"\n#####");
						}
					if(!rev.isEmpty())
						{
							System.err.println("REV is not empty! Something is BAD, BAD, BAD!");
							ArrayList<String> list = new ArrayList<String>();
							list.addAll(rev.keySet());
							Collections.sort(list);
							for (String id:list)
								{
									s = rev.get(id);
									bwSngl.write(s.toString(true));
									bwSngl.newLine();
									++snglCounter;
									System.out.println(s.getId());
								}
							System.out.println("#####\nEntries: "+rev.size()+"\n#####");
						}
					System.out.println("#########################\n#########################\nTotal of: "+snglCounter+" singletons written\n#########################\n#########################");
					
					br.close();
					bwFwd.close();
					bwRev.close();
				} 
			catch (IOException e) 
				{
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			
			
			
			
			
		}
	
	/**
	 * @param args
	 * @throws ParameterInvalidException 
	 */
	public static void main(String[] args) throws ParameterInvalidException 
		{
			String id = "$Id$";
			if(args.length==1 && args[0].equals("-v"))
				{
					System.err.println("TIGRrepairer Version: " + id);
				}
			else if(args.length==0 || (args.length==1 && args[0].equals("-h")))
				{
					TIGRrepairer.giveArgsExplanation();
				}
			else 
				{
					if(args.length>=1 && args.length <=2)
						{
							String inFile = null;
							String outDir = null;
							for(int i=0;i<args.length;++i)
								{
									if(args[i].startsWith(TIGRrepairer.INFILE_STR))
										{
											inFile=args[i].substring(TIGRrepairer.INFILE_STR.length());
										}
									if(args[i].startsWith(TIGRrepairer.OUTDIR_STR))
										outDir=args[i].substring(TIGRrepairer.OUTDIR_STR.length());
								}
							
							TIGRrepairer tr = new TIGRrepairer();
							tr.setInFile(inFile);
							tr.setOutDir(outDir);
							
							tr.readAndRePair();
						}
				}
			
		}
	
	

	

	public static void giveArgsExplanation()
		{
			System.out.println("Possible help params: -v (get version info) and -h (get this help);");
			System.out.println("Usage: java de.helmholtz_muenchen.ibis.assembly.TIGRrepairer --in=Filename [--outDir=DirectoryOut]");
			System.out.println(" --in=Filename is of course the TIGR read-file in fastq or fasta format, containing pairs of reads following the TIGR standard.");
			System.out.println(" --outDir=DirectoryOut Optional, the directory where the resulting files are to be written to. ");
			System.out.println(" if --outDir is not given, we will use the PARENT-Path of --in (Infile)");

			System.out.println("\n\n TIGR Reads (IDs) look like this:");
			System.out.println("@TIGR_C1BA001TR_926918_1103683007942\nand\n@TIGR_C1BA001TF_926922_1103671022388");
			System.out.println("The second parts are the pairing information:");
			System.out.println("C1BA001TR and C1BA001TF - R for reverse and F for forward");
			System.out.println("At least thats what we think.");
		}
	
}
