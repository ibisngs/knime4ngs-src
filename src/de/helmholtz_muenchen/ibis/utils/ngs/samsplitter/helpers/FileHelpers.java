/**
 * 
 */
package de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.helpers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.itadaki.bzip2.BZip2InputStream;
import org.itadaki.bzip2.BZip2OutputStream;


/**
 * @author jonathan.hoser
 *
 */
public class FileHelpers
	{
		public static boolean noExitRule = false;
		
		public static File checkInFile(String strInFile) throws ParameterInvalidException
			{
				return FileHelpers.checkInFile(strInFile, true);
			}
		
		public static File checkInFile(String strInFile,boolean printInfo) throws ParameterInvalidException
			{
				File inFile = new File(strInFile);
				if(!inFile.exists() || !inFile.canRead())
					{
						System.err.println("File '"+inFile.getAbsolutePath()+"' does not exist! Aborting!");
						if(FileHelpers.noExitRule)
							{
								throw new ParameterInvalidException("4: Filehelper: Given InFile does not exist.");
							}
						else
							System.exit(4);
					}
				if(printInfo)
					System.out.println("Using Infile: " + inFile.getAbsolutePath());
				return inFile;
			}
		
		//########################################################################################
		//########################################################################################
		//########################################################################################
		
		public static File checkOutDir(String strOutDir) throws ParameterInvalidException
			{
				return FileHelpers.checkOutDir(strOutDir, true);
			}
		
		public static File checkOutDir(String strOutDir,boolean printInfo) throws ParameterInvalidException
			{
				File outDir = new File(strOutDir);
				if(outDir.exists() && outDir.isDirectory())
					{
						if(printInfo)
							System.out.println("Using OutDir: " + outDir.getAbsolutePath());
						return outDir;
					}
				else
					{
						if(!outDir.exists())
							{
								if(printInfo)
									System.out.println("Creating OutDir: " + outDir.getAbsolutePath());
								outDir.mkdirs();
								if(printInfo)
									System.out.println("Using OutDir: " + outDir.getAbsolutePath());
								return outDir;
							}
						if(!outDir.isDirectory())
							{
								System.err.println("Given OutDir is not a directory: " + outDir.getAbsolutePath() + "! ABORTING!");
								if(FileHelpers.noExitRule)
									{
										throw new ParameterInvalidException("3: Filehelper: Given OutDir is not a directory.");
									}
								else
									System.exit(3);
							}
					}
				
				return outDir;
			}
		
		/**
		 * Check for a given File.
		 * If it is an existing Directory, create a new File according to strOutFileStub eg. [Dir]/[outFileStub]
		 * if it is a File, check for existance. And use it.
		 * @param strOutDir
		 * @param strOutFileStub
		 * @param printInfo
		 * @param failIfExists
		 * @return
		 * @throws ParameterInvalidException 
		 */
		public static File checkOutDirFileStub(String strOutDir, String strOutFileStub, boolean printInfo, boolean failIfExists) throws ParameterInvalidException
			{
				File outDir = new File(strOutDir);
				if(outDir.exists() && outDir.isDirectory())
					{
						File newFile = FileHelpers.checkOutFile(outDir.getAbsolutePath()+System.getProperty("file.separator")+strOutFileStub,true);
						
						return newFile;
					}
				else
					{
						if(!outDir.exists())
							{
								
								if(printInfo)
									System.out.println("Using OutFile: " + outDir.getAbsolutePath());
								return outDir;
							}
						if(!outDir.isDirectory())
							{
								System.err.println("Given OutFile exists: " + outDir.getAbsolutePath() + "!!");
								if(failIfExists)
									{
										System.err.println("Aborting -- we do NOT overwrite!");
										{
											//System.exit(1);
											if(FileHelpers.noExitRule)
												{
													throw new ParameterInvalidException("1: Filehelper: Not overwriting an outfile.");
												}
											else
												System.exit(1);
												
										}
									}
								else
									{
										System.err.println("Overwriting the existing file ....");
									}
							}
					}
				
				return outDir;
			}
		
		//########################################################################################
		//########################################################################################
		//########################################################################################
		
		public static File checkInDir(String strInDir) throws ParameterInvalidException
			{
				return FileHelpers.checkInDir(strInDir, true);
			}
		
		public static File checkInDir(String strInDir,boolean printInfo) throws ParameterInvalidException
			{
				File inDir = new File(strInDir);
				if(inDir.exists() && inDir.isDirectory())
					{
						if(printInfo)
							System.out.println("Using InDir: " + inDir.getAbsolutePath());
						return inDir;
					}
				else
					{
						if(!inDir.exists())
							{
								if(printInfo)
									System.out.println("Given InDir " + inDir.getAbsolutePath() + " does not exist! ABORTING!");
								return null;
							}
						if(!inDir.isDirectory())
							{
								System.err.println("Given InDir is not a directory: " + inDir.getAbsolutePath() + "! ABORTING!");
								
								{
									//System.exit(1);
									if(FileHelpers.noExitRule)
										{
											throw new ParameterInvalidException("2: Filehelper: Given InDir is not a directory.");
										}
									else
										System.exit(2);
										
								}
							}
					}
				
				return inDir;
			}
		
		
		
		//########################################################################################
		//########################################################################################
		//########################################################################################
		
		public static File checkOutFile(String strOutFile) throws ParameterInvalidException
			{
				return FileHelpers.checkOutFile(strOutFile,false);
			}
		
		public static File checkOutFile(String strOutFile,boolean failOnExist) throws ParameterInvalidException
			{
				return FileHelpers.checkOutFile(strOutFile, failOnExist, true);
			}

		public static File checkOutFile(String strOutFile,boolean failOnExist,boolean printInfo) throws ParameterInvalidException
			{
				File outFile = new File(strOutFile);
				if(outFile.exists())
					{
						System.err.println("File '"+outFile.getAbsolutePath()+"' exist! Overwriting!");
						if(failOnExist)
						{
							//System.exit(1);
							if(FileHelpers.noExitRule)
								{
									throw new ParameterInvalidException("1: Filehelper: Outfile exists, failing.");
								}
							else
								System.exit(1);
								
						}
					}
				if(printInfo)
					System.out.println("Using Outfile: " + outFile.getAbsolutePath());
				
				return outFile;
			}

		/**
		 * Get a BufferedReader GZ/BZ2/normal depending on the filename-ending;
		 * 
		 * @param forwardReadFile
		 * @return the WrappedBufferedReader!
		 * @throws IOException 
		 * @throws FileNotFoundException 
		 */
		public static BufferedReader getBufferedReaderOnFileEnding(File file) throws FileNotFoundException, IOException
			{
				BufferedReader brOut;
				
				if ( file.getName().toLowerCase().endsWith(".gz") )
					{
						brOut = new BufferedReader(new InputStreamReader(
								new GZIPInputStream(new FileInputStream(file))));
					}
				else if ( file.getName().toLowerCase().endsWith(".bz2") )
					{
						brOut = new BufferedReader(new InputStreamReader(
								new BZip2InputStream(new FileInputStream(file), false)));
					}
				else
					{
						brOut = new BufferedReader(new FileReader(file));
					}
				return brOut;
			}

		public static BufferedWriter getBufferedWriterOnFile(File file,
				String format) throws FileNotFoundException, IOException
			{
				BufferedWriter bw;
				
				if(format.toLowerCase().equals("auto"))
					{
						format = file.getName();
					}
				
				if ( format.toLowerCase().endsWith("gz") )
					{
						bw = new BufferedWriter(new OutputStreamWriter(
								new GZIPOutputStream(new FileOutputStream(file))));
					}
//				else if(format.toLowerCase().endsWith("bz2"))
//					{
//						bw = new BufferedWriter(new OutputStreamWriter(
//								new BZip2OutputStream(new FileOutputStream(file))));
//					}
				else
					{
						bw = new BufferedWriter(new FileWriter(file));
					}
				return bw;
				
			}
		
		
	}
