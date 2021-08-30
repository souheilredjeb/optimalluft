//@Author Souheil Rejeb
package com.optimalluft.simplexe;
public class Simplexe_Phase_II {
	
	private int m;
	private int n;
	private double[] b_vector=new double[m];
	private double[] c_vector = new double[n];
	private int[] base= new int[m];
	private double[] c_b=new double[m];
	private double[][] initial_matrix=new double[m][n];
	private Integer pivot_column_index;
	private Integer pivot_line_index;
	private double[][] augmented_matrix=new double[m+2][m+n+1];	
	private double[] pivot_column=new double[m];
	private double[] bk_column=new double[m];

	public Simplexe_Phase_II()
	{
		this.m=3;
		this.n=3;
		double[] b = {340.0,2400.0,560.0};
		this.b_vector=b;
		double[] c= {-1100.0,-1400.0,-1500.0};
		this.c_vector=c;
		double[][] matrix= {{1.0,1.0,1.0},{2.0,3.0,1.0},{1.0,2.0,3.0}};
		this.initial_matrix=matrix;
		
		int[] base=new int[this.getM()];
		for(int i=0;i<this.getM();i++)
		{
			base[i]=i+this.getN();
		}
		this.setBase(base);
		
		double[] c_b=new double[m];
		for(int i=0;i<this.getM();i++)
		{
			if(this.getBase()[i]<this.getN())
			{
				c_b[i]=this.getC_vector()[this.getBase()[i]];
			}
		}
		this.setC_b(c_b);
		
		double[][] augmented_matrix=new double[m+2][m+n+1];
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
				augmented_matrix[i][j]=1.0;
			}
		}
		for(int i=0;i<this.getM();i++)
		{
			for(int j=this.getM()+this.getN();j<this.getM()+this.getN()+1;j++)
			{
				augmented_matrix[i][j]=this.getB_vector()[i];
			}
		}
		for(int i=this.getM();i<this.getM()+1;i++)
		{
			for(int j=0;j<this.getN();j++)
			{
				augmented_matrix[i][j]=this.getC_vector()[j];
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
		this.setAugmented_matrix(augmented_matrix);
	}
	
	public Integer compute_pivot_column_index()
	{
		double[] z_c = new double[this.getM()+this.getN()];
		double max_z_c =Double.MIN_VALUE;
		int pivot_column_index=Integer.MIN_VALUE;
		for(int i=this.getM()+1;i<this.getM()+2;i++)
		{
			for(int j=0;j<this.getM()+this.getN();j++)
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
			bk_column[i]=augmented_matrix[i][this.getM()+this.getN()];
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
		for(int i=0;i<this.getM();i++)
		{
			if(i==this.getPivot_line_index())
			{
				this.getBase()[i]=this.getPivot_column_index();		
			}
		}		
	}	
	
	public void update_c_b()
	{
		double[] c_b=new double[m];
		for(int i=0;i<this.getM();i++)
		{
			if(this.getBase()[i]<this.getN())
			{
				c_b[i]=this.getC_vector()[this.getBase()[i]];
			}
		}
		this.setC_b(c_b);
	}
	
	public void apply_gauss_pivot()
	{
		double[][] augmented_matrix=new double[m+2][m+n+1];
		
		int i0=this.getPivot_line_index();
		int j0=this.getPivot_column_index();
		
		for(int i=0;i<=this.getM();i++)
		{
			if(i==i0)
			{
				for(int j=0;j<this.getM()+this.getN()+1;j++)
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
				for(int j=0;j<this.getM()+this.getN()+1;j++)
				{
					augmented_matrix[i][j]=this.getAugmented_matrix()[i][j]-this.getAugmented_matrix()[i][j0]/this.getAugmented_matrix()[i0][j0]*this.getAugmented_matrix()[i0][j];
				}
			}
		}
		for(int i=this.getM();i<this.getM()+1;i++)
		{
			for(int j=0;j<this.getN();j++)
			{
				augmented_matrix[i][j]=this.getC_vector()[j];
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
		this.setAugmented_matrix(augmented_matrix);
	}
	
	public void apply_simplexe_phase_II()
	{
		if(this.compute_pivot_column_index() >=0)
		{	
			this.build_pivot_column();
			this.build_bk_column();
			this.compute_pivot_line_index();
			this.update_base();
			this.update_c_b();
			this.apply_gauss_pivot();
			this.apply_simplexe_phase_II();
		}
	}
	
	public void display_double_vector(double[] double_vector)
	{
		for(int i=0; i< double_vector.length;i++)
		{
			System.out.println(double_vector[i]);
		}
	}
	public void display_int_vector(int[] int_vector)
	{
		for(int i=0;i<int_vector.length;i++)
		{
			System.out.println(int_vector[i]);
		}
	}	
	public void display_double_matrix(double[][] double_matrix)
	{
		for(int i=0;i<double_matrix.length;i++)
		{
			for(int j=0; j<double_matrix[0].length;j++)
			{
				System.out.print("##"); System.out.print(double_matrix[i][j]);
			}
			System.out.println();
		}
	}	
	public void initialize_base()
	{
		int[] base=this.getBase();
		for(int i=0;i<this.getM();i++)
		{
			base[i]=i+this.getN();
		}
		this.setBase(base);
	}
	
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Simplexe_Phase_II s1= new Simplexe_Phase_II();
		s1.apply_simplexe_phase_II();
		s1.display_double_matrix(s1.getAugmented_matrix());
		System.out.println("7777777777777777777777777777777777777777777777777");
		s1.display_int_vector(s1.getBase());
		
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

	public double[] getB_vector() {
		return b_vector;
	}

	public void setB_vector(double[] b_vector) {
		this.b_vector = b_vector;
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

	public Integer getPivot_column_index() {
		return pivot_column_index;
	}

	public void setPivot_column_index(Integer pivot_column_index) {
		this.pivot_column_index = pivot_column_index;
	}

	public Integer getPivot_line_index() {
		return pivot_line_index;
	}

	public void setPivot_line_index(Integer pivot_line_index) {
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
	
}
