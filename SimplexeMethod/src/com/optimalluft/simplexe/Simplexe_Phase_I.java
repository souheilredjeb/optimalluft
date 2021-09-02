//@Author Souheil Rejeb
package com.optimalluft.simplexe;


import java.util.ArrayList;

public class Simplexe_Phase_I
{
	private int m;
	private int n;
	private double[] b_vector_init=new double[m+1];
	private double[] b_vector=new double[m+1];
	private int negative_multiplicity;
	private ArrayList<Double> negative_values;
	private double[] c_vector=new double[this.m+this.n+this.negative_multiplicity];
	

	private double[] c_vector_objective=new double[this.n+this.m];
	private int[] base= new int[m];
	private double[] c_b=new double[m];
	private double[][] initial_matrix=new double[m][n];
	private int pivot_column_index;
	private int pivot_line_index;
	private int pivot_column_index_II;
	private int pivot_line_index_II;
	private double[][] augmented_matrix=new double[m+2][m+n+negative_multiplicity+1];
	private double[][] augmented_matrix_II=new double[m+2][m+n+1];


	private double[] pivot_column=new double[m];
	private double[] bk_column=new double[m];
	
	
	
	public Simplexe_Phase_I() 
	{
		// TODO Auto-generated constructor stub
		this.m=2;
		this.n=3;
		double[] b= {-6.0,-10.0};
		this.setB_vector_init(b);
		double[][] i_m= {{1.0,1.0,2.0},{1.0,2.0,1.0}};
		this.setInitial_matrix(i_m);
		
	}
	
	public void compute_negative_values()
	{
		for(int i=0;i<b_vector_init.length;i++)
		{
			if(this.b_vector_init[i]<0)
			{
				negative_multiplicity++;	
			}
		}
		this.setNegative_multiplicity(negative_multiplicity);
	}
	
	public void build_negative_values()
	{
		ArrayList<Double> n_v_=new ArrayList<Double>();
		for(int i=0;i<b_vector_init.length;i++)
		{
			if(this.b_vector_init[i]<0)
			{
				n_v_.add(b_vector_init[i]);	
			}
		}
		this.setNegative_values(n_v_);
	}
	
	public void build_c_vector()
	{
		double[] c=new double[this.m+this.n+this.negative_multiplicity];
		for(int i=0;i<this.m+this.n;i++)
		{
			c[i]=0.0;
		}
		for(int i=this.m+this.n;i<this.m+this.n+this.getNegative_multiplicity();i++)
		{
			c[i]=1.0;
		}
		this.setC_vector(c);
	}
	
	public void build_b_vector()
	{
		double[] b =new double[this.getM()];
		for(int i=0;i<this.getM();i++)
		{
			if(this.getB_vector_init()[i]<0)
			{
				b[i]=-this.getB_vector_init()[i];
			}
			else
			{
				b[i]=this.getB_vector_init()[i];
			}
		}
		this.setB_vector(b);
	}
	
	public void initialize_base()
	{
		int[] b= new int[m];
		for(int i=0;i<this.m;i++)
		{
			b[i]=this.getM()+this.getN()+i;
		}
		this.setBase(b);
	}
	
	public void initialize_c_b()
	{
		double[] c= new double[m];
		for(int i=0;i<this.getM();i++)
		{
			c[i]=1.0;
		}
		this.setC_b(c);
	}
	
	public void build_augmented_matrix()
	{
		double[][] augmented_matrix=new double[this.getM()+2][this.getM()+this.getN()+this.getNegative_multiplicity()+1];
		for(int i=0;i<this.getM();i++)
		{
			for(int j=0;j<this.getN();j++)
			{
				augmented_matrix[i][j]=this.getInitial_matrix()[i][j];
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			for(int j=this.getN();j<this.getM()+this.getN();j++)
			{
				if(j==this.getN()+i)
				{
					augmented_matrix[i][j]=-1.0;
				}
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			for(int j=this.getM()+this.getN();j<this.getM()+this.getN()+this.getNegative_multiplicity();j++)
			{
				if(j==this.getM()+this.getN()+i)
				{
					augmented_matrix[i][j]=1.0;
				}
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			for(int j=this.getM()+this.getN();j<this.getM()+this.getN()+this.getNegative_multiplicity();j++)
			{
				if(j==this.getM()+this.getN()+i)
				{
					augmented_matrix[i][j]=1.0;
				}
			}
		}
		for(int i=this.getM();i<this.getM()+1;i++)
		{
			for(int j=0;j<this.m+this.n+this.getNegative_multiplicity();j++)
			{
				augmented_matrix[i][j]=this.getC_vector()[j];		
			}
		}
		for(int i=this.getM()+1;i<this.getM()+2;i++)
		{
			for(int j=0; j<this.getN()+this.getM()+this.getNegative_multiplicity();j++)
			{
				double zj=0.0;
				for(int k=0;k<this.getM();k++)
				{
					zj=zj+augmented_matrix[k][j]*this.getC_b()[k];
				}
				augmented_matrix[i][j]=zj-augmented_matrix[this.getM()][j];
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			for(int j=this.getM()+this.getN()+this.getNegative_multiplicity();j<this.getM()+this.getN()+this.getNegative_multiplicity()+1;j++)
			{
				augmented_matrix[i][j]=this.getB_vector()[i];
			}
		}
		this.setAugmented_matrix(augmented_matrix);
		
	}
	
	
	public Integer compute_pivot_column_index()
	{
		double[] z_c = new double[this.getM()+this.getN()+this.getNegative_multiplicity()];
		double max_z_c =Double.MIN_VALUE;
		int pivot_column_index=Integer.MIN_VALUE;
		for(int i=this.getM()+1;i<this.getM()+2;i++)
		{
			for(int j=0;j<this.getM()+this.getN()+this.getNegative_multiplicity();j++)
			{
				z_c[j]=this.getAugmented_matrix()[i][j];
			}
		}
		for(int j=0;j<z_c.length;j++)
		{
			if(z_c[j] >0 && z_c[j]>max_z_c)
			{
				max_z_c=z_c[j];
				pivot_column_index=j;
			}
		}
		this.setPivot_column_index(pivot_column_index);
		return pivot_column_index;
	}
	
	public double[] build_pivot_column()
	{
		double[] pivot_column=new double[this.getM()];	
		double[][] augmented_matrix=this.getAugmented_matrix();	
		for(int i=0; i<this.getM();i++)
		{
			pivot_column[i]=augmented_matrix[i][this.getPivot_column_index()];
		}
		this.setPivot_column(pivot_column);
		return pivot_column;
	}
	
	public double[] build_bk_column()
	{
		double[] bk_column=new double[this.getM()];
		double[][] augmented_matrix=this.getAugmented_matrix();
		for(int i=0; i<this.getM();i++)
		{
			bk_column[i]=augmented_matrix[i][this.getM()+this.getN()+this.getNegative_multiplicity()];
		}
		this.setBk_column(bk_column);
		return bk_column;
	}
	
	public Integer compute_pivot_line_index()
	{
		int pivot_line_index=Integer.MIN_VALUE;
		Double[] bk_div_aik=new Double[this.getM()];
		for(int i=0;i<this.getM();i++)
		{
			if(build_bk_column()[i]>0 && this.build_pivot_column()[i]>0)
			{
				bk_div_aik[i]=this.build_bk_column()[i]/this.build_pivot_column()[i];
			}
		}
		
		for(int i=0;i<this.getM();i++)
		{
			if(bk_div_aik[i] !=null && bk_div_aik[i]>0 )
			{
				pivot_line_index=i;
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			if(bk_div_aik[i] !=null && bk_div_aik[i]<bk_div_aik[pivot_line_index])
			{
				pivot_line_index=i;
			}
		}
		this.setPivot_line_index(pivot_line_index);
		return pivot_line_index;
	}
	
	public void update_base()
	{
		int[] b=new int[this.getM()];
		for(int i=0;i<this.getM();i++)
		{
			if(i==this.getPivot_line_index())
			{
				b[i]=this.getPivot_column_index();
			}
			else
			{
				b[i]=this.getBase()[i];
			}
		}
		this.setBase(b);
	}	
	
	public void update_c_b()
	{
		double[] c_b=new double[m];
		for(int i=0;i<this.getM();i++)
		{
			if(this.getBase()[i]<this.getN()+this.getM()+this.getNegative_multiplicity())
			{
				c_b[i]=this.getC_vector()[this.getBase()[i]];
			}
		}
		this.setC_b(c_b);
	}
	
	public void apply_gauss_pivot()
	{
		double[][] augmented_matrix=new double[m+2][m+n+this.getNegative_multiplicity()+1];
		
		int i0=this.getPivot_line_index();
		int j0=this.getPivot_column_index();
		
		for(int i=0;i<=this.getM();i++)
		{
			if(i==i0)
			{
				for(int j=0;j<this.getM()+this.getN()+this.getNegative_multiplicity()+1;j++)
				{
					if( this.getAugmented_matrix()[i0][j0]!=0.0)
					{
						augmented_matrix[i][j]=this.getAugmented_matrix()[i][j]/this.getAugmented_matrix()[i0][j0];		
					}		
				}
				break;
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			if(i!=i0)
			{
				for(int j=0;j<this.getM()+this.getN()+this.getNegative_multiplicity()+1;j++)
				{
					augmented_matrix[i][j]=this.getAugmented_matrix()[i][j]-this.getAugmented_matrix()[i][j0]/this.getAugmented_matrix()[i0][j0]*this.getAugmented_matrix()[i0][j];
				}
			}
		}
		for(int i=this.getM();i<this.getM()+1;i++)
		{
			for(int j=0;j<this.m+this.n+this.getNegative_multiplicity();j++)
			{
				augmented_matrix[i][j]=this.getC_vector()[j];		
			}
		}
		for(int i=this.getM()+1;i<this.getM()+2;i++)
		{
			for(int j=0; j<this.getN()+this.getM()+this.getNegative_multiplicity();j++)
			{
				double zj=0.0;
				for(int k=0;k<this.getM();k++)
				{
					zj=zj+augmented_matrix[k][j]*this.getC_b()[k];
				}
				augmented_matrix[i][j]=zj-augmented_matrix[this.getM()][j];
			}
		}
		this.setAugmented_matrix(augmented_matrix);
	}
	
	public void apply_simplexe_phase()
	{
		if(this.compute_pivot_column_index() >=0)
		{	
			this.build_pivot_column();
			this.build_bk_column();
			this.compute_pivot_line_index();
			this.update_base();
			this.update_c_b();
			this.apply_gauss_pivot();
			this.apply_simplexe_phase();
		}
	}
	
	public void update_base_I()
	{
		int[] b=new int[this.getM()];
		for(int i=0;i<this.getM();i++)
		{
			if(i==this.getPivot_line_index())
			{
				b[i]=this.getPivot_column_index();
			}
			else
			{
				b[i]=this.getBase()[i];
			}
		}
		this.setBase(b);
	}	

	public void build_c_vector_objective()
	{
		double[] c= {1.0,1.5,3.0};
		double[] c1=new double[this.getN()+this.getM()+this.negative_multiplicity];
		for(int i=0;i<this.getN()+this.getM();i++)
		{
			if( i<this.getN())
			{
				c1[i]=c[i];
			}
			else if(i>=this.getN()+this.getM())
			{
				;
			}
			else
			{
				c1[i]=0.0;
			}
		}
		this.setC_vector_objective(c1);		
	}
	
	public void reevaluate_c_b()
	{
		double[] c_b=new double[m];
		for(int i=0;i<this.getM();i++)
		{
			if(this.getBase()[i]<this.getN()+this.getM()+this.getNegative_multiplicity())
			{
				System.out.println(this.getBase()[i]);
				c_b[i]=this.getC_vector_objective()[this.getBase()[i]];
			}
		}
		this.setC_b(c_b);
	}
	
	public void reevaluate_augmented_matrix()
	{
		double[][] augmented_matrix=new double[m+2][m+n+1];
		for(int i=this.getM();i<this.getM()+1;i++)
		{
			for(int j=0;j<this.getN()+this.getM();j++)
			{
				augmented_matrix[i][j]=this.getC_vector_objective()[j];
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			for(int j=0;j<this.getN()+this.getM();j++)
			{
				augmented_matrix[i][j]=this.getAugmented_matrix()[i][j];
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			for(int j=this.getN()+this.getM();j<this.getN()+this.getM()+1;j++)
			{
				augmented_matrix[i][j]=this.getAugmented_matrix()[i][this.getN()+this.getM()+this.getNegative_multiplicity()];
			}
		}	
		for(int i=this.getM()+1;i<this.getM()+2;i++)
		{	
			for(int j=0; j<this.getN()+this.getM();j++)
			{
				double zj=0.0;
				for(int k=0;k<this.getM();k++)
				{
					zj=zj+augmented_matrix[k][j]*this.getC_b()[k];
				}
				augmented_matrix[i][j]=zj-augmented_matrix[this.getM()][j];
			}
		}
		this.setAugmented_matrix_II(augmented_matrix);
	}
	
	public Integer compute_pivot_column_index_II()
	{
		double[] z_c = new double[this.getM()+this.getN()];
		double max_z_c =Double.MIN_VALUE;
		int pivot_column_index=Integer.MIN_VALUE;
		for(int i=this.getM()+1;i<this.getM()+2;i++)
		{
			for(int j=0;j<this.getM()+this.getN();j++)
			{
				z_c[j]=this.getAugmented_matrix_II()[i][j];
			}
		}
		for(int j=0;j<z_c.length;j++)
		{
			if(z_c[j] >0 && z_c[j]>max_z_c)
			{
				max_z_c=z_c[j];
				pivot_column_index=j;
			}
		}
		this.setPivot_column_index_II(pivot_column_index);
		return pivot_column_index;
	}
	
	public double[] build_pivot_column_II()
	{
		double[] pivot_column=new double[this.getM()];	
		double[][] augmented_matrix=this.getAugmented_matrix_II();
		
		for(int i=0; i<this.getM();i++)
		{
			pivot_column[i]=augmented_matrix[i][this.getPivot_column_index_II()];
		}
		this.setPivot_column(pivot_column);
		return pivot_column;
	}
	
	public double[] build_bk_column_II()
	{
		double[] bk_column=new double[this.getM()];
		double[][] augmented_matrix=this.getAugmented_matrix_II();
		for(int i=0; i<this.getM();i++)
		{
			bk_column[i]=augmented_matrix[i][this.getM()+this.getN()];
		}
		this.setBk_column(bk_column);
		return bk_column;
	}
	
	public Integer compute_pivot_line_index_II()
	{
		int pivot_line_index=Integer.MIN_VALUE;
		Double[] bk_div_aik=new Double[this.getM()];
		for(int i=0;i<this.getM();i++)
		{
			if(this.build_bk_column_II()[i]>0 && this.build_pivot_column_II()[i]>0)
			{
				bk_div_aik[i]=this.build_bk_column_II()[i]/this.build_pivot_column_II()[i];
			}
		}
		
		for(int i=0;i<this.getM();i++)
		{
			if(bk_div_aik[i] !=null && bk_div_aik[i]>0 )
			{
				pivot_line_index=i;
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			if(bk_div_aik[i] !=null && bk_div_aik[i]<bk_div_aik[pivot_line_index])
			{
				pivot_line_index=i;
			}
		}
		this.setPivot_line_index_II(pivot_line_index);
		return pivot_line_index;
	}
	
	public void update_base_II()
	{
		for(int i=0;i<this.getM();i++)
		{
			if(i==this.getPivot_line_index_II())
			{
				this.getBase()[i]=this.getPivot_column_index_II();		
			}
		}		
	}
	
	public void update_c_b_II()
	{
		double[] c_b=new double[m];
		for(int i=0;i<this.getM();i++)
		{
			if(this.getBase()[i]<this.getN())
			{
				c_b[i]=this.getC_vector_objective()[this.getBase()[i]];
			}
		}
		this.setC_b(c_b);
	}
	
	public void apply_gauss_pivot_II()
	{
		double[][] augmented_matrix=new double[this.getM()+2][this.getM()+this.getN()+1];
		
		int i0=this.getPivot_line_index_II();
		int j0=this.getPivot_column_index_II();
		
		for(int i=0;i<this.getM();i++)
		{
			if(i==i0)
			{
				for(int j=0;j<this.getM()+this.getN()+1;j++)
				{
					if( this.getAugmented_matrix_II()[i0][j0]!=0.0)
					{
						augmented_matrix[i][j]=this.getAugmented_matrix_II()[i][j]/this.getAugmented_matrix_II()[i0][j0];		
					}		
				}
				break;
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			if(i!=i0)
			{
				for(int j=0;j<this.getM()+this.getN()+1;j++)
				{
					
					augmented_matrix[i][j]=this.getAugmented_matrix_II()[i][j]-this.getAugmented_matrix_II()[i][j0]/this.getAugmented_matrix_II()[i0][j0]*this.getAugmented_matrix_II()[i0][j];
				
				}
			}
		}
		for(int i=this.getM();i<this.getM()+1;i++)
		{
			for(int j=0;j<this.getN();j++)
			{
				augmented_matrix[i][j]=this.getC_vector_objective()[j];
			}
			for(int j=this.getN();j<this.getM()+this.getN()+1;j++)
			{
				augmented_matrix[i][j]=0.0;
			}
			
		}
		for(int i=this.getM()+1;i<this.getM()+2;i++)
		{
			for(int j=0; j<this.getN()+this.getM();j++)
			{
				double zj=0.0;
				for(int k=0;k<this.getM();k++)
				{
					zj=zj+augmented_matrix[k][j]*this.getC_b()[k];
				}
				augmented_matrix[i][j]=zj-augmented_matrix[this.getM()][j];
			}
		}
		this.setAugmented_matrix_II(augmented_matrix);
	}
	
	public void apply_simplexe_phase_II()
	{
		if(this.compute_pivot_column_index_II()>=0)
		{	
			this.build_pivot_column_II();
			this.build_bk_column_II();
			this.compute_pivot_line_index_II();
			this.update_base_II();
			this.update_c_b_II();
			this.apply_gauss_pivot_II();
			this.apply_simplexe_phase_II();
		}
	}
	
	public void display(double[][] matrix)
	{
		for(int  i=0;i<matrix.length;i++)
		{
			for(int j=0;j<matrix[0].length;j++)
			{
				System.out.print("##");
				System.out.print(matrix[i][j]);
			}
			System.out.println();
		}
	}
	
	public void display(double[] vector)
	{
		for(int i=0;i<vector.length;i++)
		{
			System.out.println(vector[i]);
		}
	}
	
	public void display (int[] vector)
	{
		for(int i=0;i<vector.length;i++)
		{
			System.out.println(vector[i]);
		}
	}
	
	
	public static void main(String[] args) 
	{
		// TODO Auto-generated method stub
		
		Simplexe_Phase_I s1=new Simplexe_Phase_I();
		s1.compute_negative_values();
		s1.build_negative_values();
		s1.build_c_vector();
		s1.build_b_vector();
		s1.initialize_base();
		s1.initialize_c_b();
		s1.build_augmented_matrix();
		s1.apply_simplexe_phase();	
		s1.build_c_vector_objective();
		s1.reevaluate_c_b();
		s1.reevaluate_augmented_matrix();
		s1.apply_simplexe_phase_II();
		s1.display(s1.getAugmented_matrix_II());
		
	}
	
	public double[][] getAugmented_matrix_II() {
		return augmented_matrix_II;
	}

	public void setAugmented_matrix_II(double[][] augmented_matrix_II) {
		this.augmented_matrix_II = augmented_matrix_II;
	}
	
	
	
	
	
	public ArrayList<Double> getNegative_values()
	{
		return negative_values;
	}

	public void setNegative_values(ArrayList<Double> negative_values) {
		this.negative_values = negative_values;
	}

	public int getM() {
		return m;
	}

	public void setM(int m) {
		this.m = m;
	}

	public int getN() {
		return n;
	}

	public void setN(int n) {
		this.n = n;
	}

	public double[] getB_vector_init() {
		return b_vector_init;
	}

	public void setB_vector_init(double[] b_vector_init) {
		this.b_vector_init = b_vector_init;
	}

	public double[] getB_vector() {
		return b_vector;
	}

	public void setB_vector(double[] b_vector) {
		this.b_vector = b_vector;
	}

	public int getNegative_multiplicity() {
		return negative_multiplicity;
	}

	public void setNegative_multiplicity(int negative_multiplicity) {
		this.negative_multiplicity = negative_multiplicity;
	}

	public double[] getC_vector() {
		return c_vector;
	}

	public void setC_vector(double[] c_vector) {
		this.c_vector = c_vector;
	}

	public int[] getBase() {
		return base;
	}

	public void setBase(int[] base) {
		this.base = base;
	}

	public double[] getC_b() {
		return c_b;
	}

	public void setC_b(double[] c_b) {
		this.c_b = c_b;
	}

	public double[][] getInitial_matrix() {
		return initial_matrix;
	}

	public void setInitial_matrix(double[][] initial_matrix) {
		this.initial_matrix = initial_matrix;
	}

	public int getPivot_column_index() {
		return pivot_column_index;
	}

	public void setPivot_column_index(int pivot_column_index) {
		this.pivot_column_index = pivot_column_index;
	}

	public int getPivot_line_index() {
		return pivot_line_index;
	}

	public void setPivot_line_index(int pivot_line_index) {
		this.pivot_line_index = pivot_line_index;
	}

	public double[][] getAugmented_matrix() {
		return augmented_matrix;
	}

	public void setAugmented_matrix(double[][] augmented_matrix) {
		this.augmented_matrix = augmented_matrix;
	}

	public double[] getPivot_column() {
		return pivot_column;
	}

	public void setPivot_column(double[] pivot_column) {
		this.pivot_column = pivot_column;
	}

	public double[] getBk_column() {
		return bk_column;
	}

	public void setBk_column(double[] bk_column) {
		this.bk_column = bk_column;
	}

	public double[] getC_vector_objective() {
		return c_vector_objective;
	}

	public void setC_vector_objective(double[] c_vector_objective) {
		this.c_vector_objective = c_vector_objective;
	}

	public int getPivot_column_index_II() {
		return pivot_column_index_II;
	}

	public void setPivot_column_index_II(int pivot_column_index_II) {
		this.pivot_column_index_II = pivot_column_index_II;
	}

	public int getPivot_line_index_II() {
		return pivot_line_index_II;
	}

	public void setPivot_line_index_II(int pivot_line_index_II) {
		this.pivot_line_index_II = pivot_line_index_II;
	}

	

}
