package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import java.sql.*;

public class HTEDBHandler {
	
	private Connection con;
	private final String HOST = "ibisdb01";
	private final String DBNAME = "ngs_HTE";
	private final String USER = "ngs_hteuser";
	private final String PW = "htepass";
	
	
	public HTEDBHandler() throws SQLException {
		try {
			Class.forName("com.mysql.jdbc.Driver").newInstance();
			this.con = DriverManager.getConnection("jdbc:mysql://"+HOST+"/"+DBNAME,USER,PW);
		} catch (InstantiationException | IllegalAccessException
				| ClassNotFoundException e) {
			e.printStackTrace();
			System.err.println("Failure while loading jdbc driver: "+e.getMessage());
		}
	}
	
	/**
	 * 
	 * @param hte_threshold
	 * @return id of inserted HTE_Workflow
	 */
	public int insertNewHTEWorkflow(int hte_threshold) {
		String query = "INSERT INTO HTE_Workflow (hte_threshold) VALUES (?)";
		PreparedStatement preStmt;
		ResultSet rs;
		int id=-1;
		try {
			preStmt = con.prepareStatement(query,Statement.RETURN_GENERATED_KEYS);
			preStmt.setInt(1,hte_threshold);
			preStmt.execute();
			rs = preStmt.getGeneratedKeys();
			if(rs.next()) {
				id = rs.getInt(1);
			}
		} catch (SQLException e) {
			System.err.println("New HTE_Workflow could not be inserted: "+e.getMessage());
		}
		return id;
	}
	
	public int insertNewHTENode(String command, String name, String host, int node_threshold, int hte_id) {
		String query = "INSERT INTO HTE_Node (klock_id, klock_string, node_name, host_name, node_threshold,hte_id) VALUES (?,?,?,?,?,?)";
		PreparedStatement preStmt;
		ResultSet rs;
		String id_string = command+name+host;
		int id = id_string.hashCode();
		try {
			preStmt = con.prepareStatement(query,Statement.RETURN_GENERATED_KEYS);
			preStmt.setInt(1,id);
			preStmt.setString(2, command);
			preStmt.setString(3, name);
			preStmt.setString(4, host);
			preStmt.setInt(5, node_threshold);
			preStmt.setInt(6, hte_id);
			preStmt.execute();
			rs = preStmt.getGeneratedKeys();
			if(rs.next()) {
				id = rs.getInt(1);
			}
		} catch (SQLException e) {
			System.err.println("New HTE_Node could not be inserted: "+e.getMessage());
		}
		return id;
	}
	
	public int getHTENodeCount(String command, String name, String host, int hte_id) {
		String id_string = command+name+host;
		int id = id_string.hashCode();
		return this.getHTENodeCount(id, hte_id);
	}
	
	public int getHTENodeCount(int node_id, int hte_id) {
		ResultSet rs = this.getHTENodeTuple(node_id, hte_id);
		try {
			if(rs.next()) {
				return rs.getInt(9);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return -1;
	}
	
	public int getHTENodeThreshold(String command, String name, String host, int hte_id) {
		String id_string = command+name+host;
		int id = id_string.hashCode();
		return this.getHTENodeCount(id, hte_id);
	}
	
	public int getHTENodeThreshold(int node_id, int hte_id) {
		ResultSet rs = this.getHTENodeTuple(node_id, hte_id);
		try {
			if(rs.next()) {
				return rs.getInt(8);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return -1;
	}
	
	public int getHTENodeId(String command, String name, String host, int hte_id) {
		String id_string = command+name+host;
		int id = id_string.hashCode();
		return this.getHTENodeId(id, hte_id);
	}
	
	public int getHTENodeId(int node_id, int hte_id) {
		ResultSet rs = this.getHTENodeTuple(node_id, hte_id);
		try {
			if(rs.next()) {
				return rs.getInt(1);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return -1;
	}
	
	private ResultSet getHTENodeTuple(int node_id, int hte_id) {
		String query = "SELECT * FROM HTE_Node WHERE klock_id = ? AND hte_id = ?";
		PreparedStatement preStmt;
		ResultSet rs = null;
		int id = node_id;
		try {
			preStmt = con.prepareStatement(query);
			preStmt.setInt(1, id);
			preStmt.setInt(2, hte_id);
			rs = preStmt.executeQuery();	
		} catch (SQLException e) {
			System.err.println("Check for klock_id "+id+" failed.");
		}
		return rs;
	}
	
	public void updateCount(int hte_id, int node_id, int count) {
		String query = "UPDATE HTE_Node SET node_count = ? WHERE hte_id = ? AND klock_id = ?";
		PreparedStatement preStmt;
		try {
			preStmt = con.prepareStatement(query);
			preStmt.setInt(1, count);
			preStmt.setInt(2, hte_id);
			preStmt.setInt(3, node_id);
			preStmt.execute();
		} catch (SQLException e) {
			System.err.println("Count could not be updated for node with klock_id "+node_id+": "+e.getMessage());
		}
	}
	
	public void updateThreshold(int hte_id, int node_id, int threshold) {
		String query = "UPDATE HTE_Node SET node_threshold = ? WHERE hte_id = ? AND klock_id = ?";
		PreparedStatement preStmt;
		try {
			preStmt = con.prepareStatement(query);
			preStmt.setInt(1, threshold);
			preStmt.setInt(2, hte_id);
			preStmt.setInt(3, node_id);
			preStmt.execute();
		} catch (SQLException e) {
			System.err.println("Count could not be updated for node with klock_id "+node_id+": "+e.getMessage());
		}
	}
	
	public void writeSuccess(int hte_id, int node_id) {
		String query = "UPDATE HTE_Node SET node_successful_finished = 1 WHERE hte_id = ? AND klock_id = ?";
		PreparedStatement preStmt;
		try {
			preStmt = con.prepareStatement(query);
			preStmt.setInt(1, hte_id);
			preStmt.setInt(2, node_id);
			preStmt.execute();
		} catch (SQLException e) {
			System.err.println("Successful could not be updated for node with klock_id "+node_id+": "+e.getMessage());
		}
	}
	
	public void writeError(int hte_id, int node_id, String errmsg) {
		String query = "INSERT INTO HTE_History (hte_id, klock_id, err_msg) VALUES (?,?,?);";
		PreparedStatement preStmt;
		try {
			preStmt = con.prepareStatement(query);
			preStmt.setInt(1, hte_id);
			preStmt.setInt(2, node_id);
			preStmt.setString(3, errmsg);
			preStmt.execute();
		} catch (SQLException e) {
			System.err.println("Error message could not be inserted for node with klock_id "+node_id+": "+e.getMessage());
		}
	}
	
	public void closeConnection() {
		try {
			this.con.close();
		} catch (SQLException e) {
			System.err.println("Connection to "+DBNAME+" couldn't be closed correctly: "+e.getMessage());
		}
	}
}
