package de.helmholtz_muenchen.ibis.utils.threads;

public class UnsuccessfulExecutionException extends Exception{

	/**
	 * 
	 */
	private static final long serialVersionUID = 5118485654421200093L;

	public UnsuccessfulExecutionException(String m){
		super(m);
	}
    public UnsuccessfulExecutionException(String message, Throwable throwable) {
        super(message, throwable);
    }
    
}
