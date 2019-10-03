
using System;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.FEM.Elements.SupportiveClasses
{
	#region

	#endregion

	/// <summary>
	/// An one-dimensional integration point generated utilizing Gauss-Legendre quadrature abscissa.
	/// </summary>
	public class GaussLegendrePoint1D
	{
		#region Properties

		/// <summary>
		/// Parametric coordinate of <see cref="GaussLegendrePoint1D"/>.
		/// </summary>
		public double Coordinate { get; set; }

		/// <summary>
		/// Weight of the <see cref="GaussLegendrePoint1D"/>.
		/// </summary>
		public double WeightFactor { get; set; }

		#endregion
	}

	/// <summary>
	/// Three-dimensional integration point.
	/// </summary>
	public class GaussLegendrePoint3D
	{
		private IMatrixView B;
		private double Ksi;
		private double Heta;

		/// <summary>
		/// Defines a <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		/// <param name="ksi">Parametric coordinate ksi of the <see cref="GaussLegendrePoint3D"/>.</param>
		/// <param name="heta">Parametric coordinate heta of the <see cref="GaussLegendrePoint3D"/>.</param>
		/// <param name="zeta">Parametric coordinate zeta of the <see cref="GaussLegendrePoint3D"/>.</param>
		/// <param name="deformationMatrix">An <see cref="IMatrixView"/> containing the deformation matrix calculated at the specific <see cref="GaussLegendrePoint3D"/>.</param>
		/// <param name="weightFactor">Weight of the <see cref="GaussLegendrePoint3D"/>.</param>
		public GaussLegendrePoint3D(double ksi, double heta, double zeta, IMatrixView deformationMatrix, double weightFactor)
		{
			this.Ksi = ksi;
			this.Heta = heta;
			this.Zeta = zeta;
			this.B = deformationMatrix;
			WeightFactor = weightFactor;
		}
		#region Constructors and Destructors

		/// <summary>
		/// Defines a <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		/// <param name="xi">Parametric coordinate ksi of the <see cref="GaussLegendrePoint3D"/>.</param>
		/// <param name="eta">Parametric coordinate heta of the <see cref="GaussLegendrePoint3D"/>.</param>
		/// <param name="zeta">Parametric coordinate zeta of the <see cref="GaussLegendrePoint3D"/>.</param>
		/// <param name="deformationMatrix">A 2D array containing the deformation matrix calculated at the specific <see cref="GaussLegendrePoint3D"/>.</param>
		/// <param name="weightFactor">Weight of the <see cref="GaussLegendrePoint3D"/>.</param>
		public GaussLegendrePoint3D(
			double xi, double eta, double zeta, double[,] deformationMatrix, double weightFactor)
		{
			this.Xi = xi;
			this.Eta = eta;
			this.Zeta = zeta;
			this.DeformationMatrix = deformationMatrix;
			this.WeightFactor = weightFactor;
		}

		#endregion

		#region Properties

		/// <summary>
		/// Returns the deformation matrix of the <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		public double[,] DeformationMatrix { get; private set; }

		/// <summary>
		/// Returns parametric coordinate Heta of the <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		public double Eta { get; private set; }

		/// <summary>
		/// Returns the weight factor of the <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		public double WeightFactor { get; private set; }

		/// <summary>
		/// Returns the parametric coordinate Ksi of the <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		public double Xi { get; private set; }

		/// <summary>
		/// Returns the parametric coordinate Zeta of the <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		public double Zeta { get; private set; }

		#endregion
	}

	/// <summary>
	/// This class defines the basics of Gauss-Legendre quadrature.
	/// </summary>
	public class GaussQuadrature
	{
		#region Constants and Fields

		private static readonly GaussLegendrePoint1D GaussLegendrePoint1 = new GaussLegendrePoint1D
		{
			Coordinate = 0.0,
			WeightFactor = 2.0
		};

		private static readonly GaussLegendrePoint1D GaussLegendrePoint2A = new GaussLegendrePoint1D
		{
			Coordinate = -0.5773502691896,
			WeightFactor = 1.0
		};

		private static readonly GaussLegendrePoint1D GaussLegendrePoint2B = new GaussLegendrePoint1D
		{
			Coordinate = 0.5773502691896,
			WeightFactor = 1.0
		};

		private static readonly GaussLegendrePoint1D GaussLegendrePoint3A = new GaussLegendrePoint1D
		{
			Coordinate = -0.7745966692415,
			WeightFactor = 0.5555555555556
		};

		private static readonly GaussLegendrePoint1D GaussLegendrePoint3B = new GaussLegendrePoint1D
		{
			Coordinate = 0.0,
			WeightFactor = 0.8888888888889
		};

		private static readonly GaussLegendrePoint1D GaussLegendrePoint3C = new GaussLegendrePoint1D
		{
			Coordinate = 0.7745966692415,
			WeightFactor = 0.5555555555556
		};

		private static readonly GaussLegendrePoint1D GaussLegendrePoint4A = new GaussLegendrePoint1D
		{
			Coordinate = -0.86113631159416,
			WeightFactor = 0.3478548451375
		};

		private static readonly GaussLegendrePoint1D GaussLegendrePoint4B = new GaussLegendrePoint1D
		{
			Coordinate = -0.3399810435849,
			WeightFactor = 0.6521451548625
		};

		private static readonly GaussLegendrePoint1D GaussLegendrePoint4C = new GaussLegendrePoint1D
		{
			Coordinate = 0.3399810435849,
			WeightFactor = 0.6521451548625
		};

		private static readonly GaussLegendrePoint1D GaussLegendrePoint4D = new GaussLegendrePoint1D
		{
			Coordinate = 0.86113631159416,
			WeightFactor = 0.3478548451375
		};

		#endregion

		#region Public Methods

		/// <summary>
		/// Returns the one dimensional integration points for a given degree.
		/// </summary>
		/// <param name="integrationDegree">Degree of integrations. The value is equal to the number of the integration points returned.</param>
		/// <returns>An array of one dimensional integration points.</returns>
		public static GaussLegendrePoint1D[] GetGaussLegendrePoints(int integrationDegree)
		{
			if (integrationDegree < 1)
			{
				throw new InvalidOperationException("Integration Degree must be greater or equal to 1. ");
			}

			switch (integrationDegree)
			{
				case 1:
					return new[] { GaussLegendrePoint1 };
				case 2:
					return new[] { GaussLegendrePoint2A, GaussLegendrePoint2B };
				case 3:
					return new[] { GaussLegendrePoint3A, GaussLegendrePoint3B, GaussLegendrePoint3C };
				case 4:
					return new[]
						{
						   GaussLegendrePoint4A, GaussLegendrePoint4B, GaussLegendrePoint4C, GaussLegendrePoint4D
						};
				default:
					throw new NotImplementedException("Integration Degree higher than 4 is not implemented yet. ");
			}
		}

		#endregion
	}
}
