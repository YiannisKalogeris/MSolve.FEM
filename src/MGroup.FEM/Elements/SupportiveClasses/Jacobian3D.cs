using System;

namespace MGroup.FEM.Elements.SupportiveClasses
{
	/// <summary>
	/// Calculates the jacobian for a given integration point.
	/// </summary>
	public class Jacobian3D
	{
		#region Constants and Fields

		private const int Dimensions = 3;

		#endregion

		#region Constructors and Destructors

		/// <summary>
		/// Defines a <see cref="Jacobian3D"/> object.
		/// </summary>
		/// <param name="nodeCoordinates">The nodal coordinates of the element.</param>
		/// <param name="shapeFunctionDerivatives">The shape functions natural derivatives of a given integration point.</param>
		public Jacobian3D(double[,] nodeCoordinates, ShapeFunctionNaturalDerivatives3D[] shapeFunctionDerivatives)
		{
			this.Matrix = CalculateJacobianMatrix(nodeCoordinates, shapeFunctionDerivatives);
			this.Determinant = this.CalculateJacobianDeterminant();
		}

		#endregion

		#region Properties

		/// <summary>
		/// Returns the jacobian matrix determinant.
		/// </summary>
		public double Determinant { get; private set; }

		/// <summary>
		/// Returns the jacobian matrix.
		/// </summary>
		public double[,] Matrix { get; private set; }

		#endregion

		#region Public Methods

		/// <summary>
		/// Calculates the inverse of teh jacobian matrix.
		/// </summary>
		/// <returns>A 2D double array containing the jacobian inverse values.</returns>
		public double[,] CalculateJacobianInverse()
		{
			double[,] jacobianInverse = new double[Dimensions, Dimensions];
			double determinantInverse = 1.0 / this.Determinant;

			jacobianInverse[0, 0] = ((this.Matrix[1, 1] * this.Matrix[2, 2]) - (this.Matrix[2, 1] * this.Matrix[1, 2])) *
									determinantInverse;
			jacobianInverse[0, 1] = ((this.Matrix[2, 1] * this.Matrix[0, 2]) - (this.Matrix[0, 1] * this.Matrix[2, 2])) *
									determinantInverse;
			jacobianInverse[0, 2] = ((this.Matrix[0, 1] * this.Matrix[1, 2]) - (this.Matrix[1, 1] * this.Matrix[0, 2])) *
									determinantInverse;
			jacobianInverse[1, 0] = ((this.Matrix[2, 0] * this.Matrix[1, 2]) - (this.Matrix[1, 0] * this.Matrix[2, 2])) *
									determinantInverse;
			jacobianInverse[1, 1] = ((this.Matrix[0, 0] * this.Matrix[2, 2]) - (this.Matrix[2, 0] * this.Matrix[0, 2])) *
									determinantInverse;
			jacobianInverse[1, 2] = ((this.Matrix[1, 0] * this.Matrix[0, 2]) - (this.Matrix[0, 0] * this.Matrix[1, 2])) *
									determinantInverse;
			jacobianInverse[2, 0] = ((this.Matrix[1, 0] * this.Matrix[2, 1]) - (this.Matrix[2, 0] * this.Matrix[1, 1])) *
									determinantInverse;
			jacobianInverse[2, 1] = ((this.Matrix[2, 0] * this.Matrix[0, 1]) - (this.Matrix[2, 1] * this.Matrix[0, 0])) *
									determinantInverse;
			jacobianInverse[2, 2] = ((this.Matrix[0, 0] * this.Matrix[1, 1]) - (this.Matrix[1, 0] * this.Matrix[0, 1])) *
									determinantInverse;
			return jacobianInverse;
		}

		#endregion

		#region Methods

		private static double[,] CalculateJacobianMatrix(
			double[,] nodeCoordinates, ShapeFunctionNaturalDerivatives3D[] shapeFunctionDerivatives)
		{
			double[,] jacobianMatrix = new double[Dimensions, Dimensions];
			for (int i = 0; i < nodeCoordinates.GetLength(0); i++)
			{
				for (int j = 0; j < Dimensions; j++)
				{
					jacobianMatrix[0, j] += nodeCoordinates[i, j] * shapeFunctionDerivatives[i].Xi;
					jacobianMatrix[1, j] += nodeCoordinates[i, j] * shapeFunctionDerivatives[i].Eta;
					jacobianMatrix[2, j] += nodeCoordinates[i, j] * shapeFunctionDerivatives[i].Zeta;
				}
			}

			return jacobianMatrix;
		}

		private double CalculateJacobianDeterminant()
		{
			double det1 = this.Matrix[0, 0] *
						  ((this.Matrix[1, 1] * this.Matrix[2, 2]) - (this.Matrix[2, 1] * this.Matrix[1, 2]));
			double det2 = this.Matrix[0, 1] *
						  ((this.Matrix[1, 0] * this.Matrix[2, 2]) - (this.Matrix[2, 0] * this.Matrix[1, 2]));
			double det3 = this.Matrix[0, 2] *
						  ((this.Matrix[1, 0] * this.Matrix[2, 1]) - (this.Matrix[2, 0] * this.Matrix[1, 1]));

			double jacobianDeterminant = det1 - det2 + det3;

			if (jacobianDeterminant < 0)
			{
				throw new InvalidOperationException("The Jacobian Determinant is negative.");
			}

			return jacobianDeterminant;
		}

		#endregion
	}
}
