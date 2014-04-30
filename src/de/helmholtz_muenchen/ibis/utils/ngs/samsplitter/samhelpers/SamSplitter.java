/**
 * 
 */
package de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.samhelpers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.helpers.FileHelpers;
import de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.helpers.ParameterInvalidException;


/**
 * @author jonathan.hoser
 *
 */
public class SamSplitter
	{
		
		private static final String OUTDIR_STR="--outDir=";
		private static final String INFILE_STR="--in=";
		private static final String SPLITSIZE_STR = "--splitsize=";
		private File inFile;
		private File outFileDir;
		private File outFileStub;
		private int samSplitSize;
		public static boolean noExit = false; // If no exit is set, instead of System.exit() a Exception will be throwed
		
		public void setOutDir(String strOutDir) throws ParameterInvalidException 
			{
				if(strOutDir!=null)
					{
						File _outDir = FileHelpers.checkOutDir(strOutDir);
						
						if(_outDir.isDirectory())
							this.outFileDir = _outDir;
						else
							{
								System.err.println("We need a directory to place the files, not a file!");
								if(!noExit) 
									System.exit(5);
								else
									throw new IllegalArgumentException("We need a directory to place the files, not a file!");
							}
						
						
					}
				else
					{
						if(this.inFile == null)
							throw new ParameterInvalidException("Error 1: inFile CANNOT be empty! But it is");
						
						//grab the parent-File/Dir for OutDir.
						
						this.outFileDir = new File(this.inFile.getAbsoluteFile().getParentFile().getAbsolutePath()+System.getProperty("file.separator"));
						
					}
				System.out.println("Writing to:");
				System.out.println("DIR:  " + this.outFileDir.getAbsolutePath());
				
			}

		public void setInFile(String strInFile) throws ParameterInvalidException 
			{
				if(strInFile!=null)
					{
						File _inFile = FileHelpers.checkInFile(strInFile);
						this.inFile = _inFile;
					}
				else
					{
						System.err.println("--in=File cannot be null;");
						SamSplitter.giveArgsExplanation();
						if(!noExit) 
							System.exit(6);
						else
							throw new IllegalArgumentException("--in=File cannot be null;");
					}
			}
		
		public void setSplitSize(int splitSizeInt)
			{
				if(splitSizeInt>0)
					this.samSplitSize = splitSizeInt;
				else
					{
						System.err.println("Split-Size cannot be negative, setting to 1mio");
						this.samSplitSize = 1000000;
					}
			}
		
		
		
		public void split()
			{
				String line;
				long a = System.nanoTime();
				
				BufferedWriter bwOut = null;
				BufferedReader br;
				int splitCount=0;
				
				StringBuffer header = new StringBuffer();
				
				
				try
					{
						//@2007-12-17_201B2AAXX_celegans_2_1_435_274
						//			
						br = new BufferedReader(new FileReader(this.inFile));
						
						//1. try to identify sam header!
						
						while((line=br.readLine())!=null) //name line
							{
								//@ means header!
								if(line.charAt(0)=='@')
									{
										header.append(line);
										header.append('\n');
										continue;
									}
								else
									break;
							}
						//line is filled with first sequence-line.
						//2. start a new outfile.
						//3. write header to it.
						//4. write sam-entries, until split-size reached.
						//5. restart with 2.
						
						while( wroteToSplitSize(header, bwOut, splitCount, br, line) == true)
							{
								line=br.readLine();
								if(line != null)
									{
										++splitCount;
									}
								else
									break;
							}
						
						
						br.close();
					} 
				catch (IOException e) 
					{
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
			}
		
		private boolean wroteToSplitSize(StringBuffer header, BufferedWriter bwOut, int splitCount, BufferedReader br, String line) throws IOException
			{
				bwOut = new BufferedWriter(new FileWriter(new File(
						this.outFileDir.getAbsolutePath()+
						System.getProperty("file.separator")+
						this.inFile.getName()+
						"." + String.valueOf(splitCount)
						)));
				
				bwOut.write(header.toString());
				//last should have a line-feed.
				//now we write the line that was not a header above.
				bwOut.write(line);
				bwOut.newLine();
				int samCounter=1;
				
				while((line=br.readLine())!=null) //name line
					{
						bwOut.append(line);
						bwOut.newLine();
						++samCounter;
						
						if(samCounter>= this.samSplitSize)
							break;
					}
				bwOut.flush();
				bwOut.close();
				if(samCounter<this.samSplitSize)
					return false;
				else
					return true;
				
			}
		
		/**
		 * @param args
		 */
		public static void main(String[] args) throws ParameterInvalidException
			{
				String id = "$Id$";
				if(args.length==1 && args[0].equals("-v"))
					{
						System.err.println("SamSplitter Version: " + id);
					}
				else if(args.length==0 || (args.length==1 && args[0].equals("-h")))
					{
						SamSplitter.giveArgsExplanation();
					}
				else 
					{
						if(args.length>=2 && args.length <=3)
							{
								String inFile = null;
								String outDir = null;
								String splitSize = null;
								int splitSizeInt = 0;
								for(int i=0;i<args.length;++i)
									{
										if(args[i].startsWith(SamSplitter.INFILE_STR))
											{
												inFile=args[i].substring(SamSplitter.INFILE_STR.length());
											}
										if(args[i].startsWith(SamSplitter.OUTDIR_STR))
											outDir=args[i].substring(SamSplitter.OUTDIR_STR.length());
										if(args[i].startsWith(SamSplitter.SPLITSIZE_STR))
											{
												splitSize = args[i].substring(SamSplitter.SPLITSIZE_STR.length());
												splitSizeInt = Integer.parseInt(splitSize);
											}
									}
								
								SamSplitter ssr = new SamSplitter();
								ssr.setInFile(inFile);
								ssr.setOutDir(outDir);
								ssr.setSplitSize(splitSizeInt);
								
								ssr.split();
							}
					}
			}
		

		

		

		public static void giveArgsExplanation()
			{
				System.out.println("Possible help params: -v (get version info) and -h (get this help);");
				System.out.println("Usage: java de.helmholtz_muenchen.ibis.samhelpers.SamSplitter --in=Filename [--outDir=DirectoryOut]");
				System.out.println(" --splitsize=xx with xx being number of SAM alignments in the SAM file.");
				System.out.println(" --in=Filename is of course the SAMfile in SAM format (not BAM!!)");
				System.out.println(" --outDir=DirectoryOut Optional, the directory where the resulting files are to be written to. ");
				System.out.println(" if --outDir is not given, we will use the PARENT-Path of --in (Infile)");
		
				
			}
		
	}
