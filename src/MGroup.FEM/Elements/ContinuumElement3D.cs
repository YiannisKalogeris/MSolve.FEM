using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.FEM.Interpolation;
using MGroup.FEM.Interpolation.GaussPointExtrapolation;
using MGroup.FEM.Interpolation.Jacobians;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Materials;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.FreedomDegrees;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Interfaces;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Elements
{
	/// <summary>
	/// Represents a continuum finite element for 3D problems. Specific elements (e.g. Hexa8, Hexa20, ...) can be created using
	/// the appropriate <see cref="IIsoparametricInterpolation3D_OLD"/>, <see cref="IQuadrature3D"/> etc. strategies. 
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class ContinuumElement3D : IStructuralFiniteElement, ICell<Node>
	{
		private readonly static IDofType[] nodalDOFTypes = new IDofType[]
		{
			StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
		};

		private readonly IDofType[][] dofTypes;
		private DynamicMaterial dynamicProperties;
		private readonly IReadOnlyList<ElasticMaterial3D> materialsAtGaussPoints;


		/// <summary>
		/// Defines a three-dimensional continuum element.
		/// </summary>
		/// <param name="nodes">The element nodes.</param>
		/// <param name="interpolation"> The interpolation method of the element.</param>
		/// <param name="quadratureForStiffness">A quadrature rule for calculating the elemental stiffness matrix.</param>
		/// <param name="quadratureForMass">A quadrature rule for calculating the elemental mass matrix.</param>
		/// <param name="gaussPointExtrapolation">A rule for extrapolating Gauss point values to the element nodes.</param>
		/// <param name="materialsAtGaussPoints">a <see cref="List{T}"/> of materials for each one of the Gauss points.</param>
		/// <param name="dynamicProperties">The dynamic properties of the element.</param>
		public ContinuumElement3D(IReadOnlyList<Node> nodes, IIsoparametricInterpolation3D interpolation,
			IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
			IGaussPointExtrapolation3D gaussPointExtrapolation,
			IReadOnlyList<ElasticMaterial3D> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
		{
			this.dynamicProperties = dynamicProperties;
			this.materialsAtGaussPoints = materialsAtGaussPoints;
			this.GaussPointExtrapolation = gaussPointExtrapolation;
			this.Nodes = nodes;
			this.Interpolation = interpolation;
			this.QuadratureForConsistentMass = quadratureForMass;
			this.QuadratureForStiffness = quadratureForStiffness;

			dofTypes = new IDofType[nodes.Count][];
			for (int i = 0; i < nodes.Count; i++)
			{
				dofTypes[i] = new IDofType[]
				{
					StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
				};
			}
		}

		/// <summary>
		/// Retrieves the type of the finite element used.
		/// </summary>
		public CellType CellType => Interpolation.CellType;

		/// <summary>
		/// Defines the way that elemental degrees of freedom will be enumerated.
		/// For further info see <see cref="IElementDofEnumerator"/>.
		/// </summary>
		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		/// <summary>
		/// Dimensionality of the element.
		/// </summary>
		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		/// <summary>
		/// Returns the rule for extrapolating Gauss point values to the element nodes.
		/// </summary>
		public IGaussPointExtrapolation3D GaussPointExtrapolation { get; }

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="ContinuumElement3D"/>.</param>
		/// <returns>A <see cref="IReadOnlyList{T}"/> that contains a <see cref="IReadOnlyList{T}"/> of <see cref="IDofType"/> with degrees of freedom for each elemental <see cref="Node"/>.</returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

		/// <summary>
		/// The element ID.
		/// </summary>
		public int ID => throw new NotImplementedException(
			"Element type codes should be in a settings class. Even then it's a bad design choice");

		/// <summary>
		/// The interpolation function used for the element.
		/// </summary>
		public IIsoparametricInterpolation3D Interpolation { get; }

		/// <summary>
		/// Boolean denoting if the material of the element has been modified.
		/// </summary>
		public bool MaterialModified
		{
			get
			{
				foreach (ElasticMaterial3D material in materialsAtGaussPoints)
				{
					if (material.Modified) return true;
				}
				return false;
			}
		}

		/// <summary>
		/// A <see cref="List{T}"/> of the Nodes of the element.
		/// </summary>
		public IReadOnlyList<Node> Nodes { get; }

		/// <summary>
		/// The quadrature rule for calculating the mass matrix.
		/// </summary>
		public IQuadrature3D QuadratureForConsistentMass { get; }

		/// <summary>
		/// The quadrature rule for calculating the stiffness matrix.
		/// </summary>
		public IQuadrature3D QuadratureForStiffness { get; }

		/// <summary>
		/// Builds the element consistent mass matrix.
		/// </summary>
		/// <returns>A <see cref="Matrix"/>.</returns>
		public Matrix BuildConsistentMassMatrix()
		{
			int numberOfDofs = 3 * Nodes.Count;
			var mass = Matrix.CreateZero(numberOfDofs, numberOfDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				Matrix partial = shapeFunctionMatrix.MultiplyRight(shapeFunctionMatrix, true, false);
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				mass.AxpyIntoThis(partial, dA);
			}
			mass.ScaleIntoThis(dynamicProperties.Density);
			return mass;
		}

		/// <summary>
		/// Builds the element lumped mass matrix.
		/// </summary>
		/// <returns>A <see cref="Matrix"/>.</returns>
		public Matrix BuildLumpedMassMatrix()
		{
			int numberOfDofs = 3 * Nodes.Count;
			var lumpedMass = Matrix.CreateZero(numberOfDofs, numberOfDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			double area = 0;
			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				area += jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
			}

			double nodalMass = area * dynamicProperties.Density / Nodes.Count;
			for (int i = 0; i < numberOfDofs; i++) lumpedMass[i, i] = nodalMass;

			return lumpedMass;
		}

		/// <summary>
		/// Calculates the stiffness matrix of the element.
		/// </summary>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="ContinuumElement3D"/>.</returns>
		public IMatrix BuildStiffnessMatrix()
		{
			int numberOfDofs = 3 * Nodes.Count;
			var stiffness = Matrix.CreateZero(numberOfDofs, numberOfDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				IMatrixView constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				Matrix shapeGradientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
				Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);

				Matrix partial = deformation.ThisTransposeTimesOtherTimesThis(constitutive);
				double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
				stiffness.AxpyIntoThis(partial, dA);
			}

			return DofEnumerator.GetTransformedMatrix(stiffness);
		}

		/// <summary>
		/// Calculates the forces applies to an <see cref="ContinuumElement3D"/> due to <see cref="MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="ContinuumElement3D"/>.</param>
		/// <param name="loads">A list of <see cref="MassAccelerationLoad"/>. For more info see <seealso cref="MassAccelerationLoad"/>.</param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			int numberOfDofs = 3 * Nodes.Count;
			var accelerations = new double[numberOfDofs];
			IMatrix massMatrix = MassMatrix(element);

			foreach (var load in loads)
			{
				int index = 0;
				foreach (var nodalDOFTypes in dofTypes)
				{
					foreach (var dofType in nodalDOFTypes)
					{
						if (dofType == load.DOF) accelerations[index] += load.Amount;
						index++;
					}
				}
			}

			return massMatrix.Multiply(accelerations);
		}

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="ContinuumElement3D"/>.</param>
		/// <param name="localTotalDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForces(IElement element, double[] localTotalDisplacements, double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="ContinuumElement3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
		{
			return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
		}

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="ContinuumElement3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="Tuple{T1,T2}"/> of the stresses and strains of the element.</returns>
		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements,
			double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// Calculates the volume of the element.
		/// </summary>
		/// <returns>The volume as a <see cref="double"/>.</returns>
		public double CalculateVolume()
		{
			double volume = 0.0;
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				volume += jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
			}
			return volume;
		}

		/// <summary>
		/// Clear the material state of the element.
		/// </summary>
		public void ClearMaterialState()
		{
			foreach (var material in materialsAtGaussPoints) material.ClearState();
		}

		/// <summary>
		/// Clear any saved material stresses of the element.
		/// </summary>
		public void ClearMaterialStresses()
		{
			foreach (var material in materialsAtGaussPoints) material.ClearStresses();
		}

		/// <summary>
		/// Calculates the damping matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="ContinuumElement3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of an <see cref="ContinuumElement3D"/>.</returns>
		public IMatrix DampingMatrix(IElement element)
		{
			IMatrix damping = BuildStiffnessMatrix();
			damping.ScaleIntoThis(dynamicProperties.RayleighCoeffStiffness);
			damping.AxpyIntoThis(MassMatrix(element), dynamicProperties.RayleighCoeffMass);
			return damping;
		}

		/// <summary>
		/// Calculates the coordinates of the centroid of this element.
		/// </summary>
		public CartesianPoint FindCentroid()
			=> Interpolation.TransformNaturalToCartesian(Nodes, new NaturalPoint(0.0, 0.0, 0.0));

		/// <summary>
		/// Calculates the mass matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="ContinuumElement3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="ContinuumElement3D"/>.</returns>
		public IMatrix MassMatrix(IElement element)
		{
			return BuildLumpedMassMatrix();
		}

		/// <summary>
		/// Resets any saved material states of the element to its initial state.
		/// </summary>
		public void ResetMaterialModified()
		{
			foreach (var material in materialsAtGaussPoints) material.ResetModified();
		}

		/// <summary>
		/// Save the current material state of the element.
		/// </summary>
		public void SaveMaterialState()
		{
			foreach (var m in materialsAtGaussPoints) m.SaveState();
		}

		/// <summary>
		/// Calculates the stiffness matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="ContinuumElement3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="ContinuumElement3D"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element) => DofEnumerator.GetTransformedMatrix(BuildStiffnessMatrix());


		/// <summary>
		/// Updates the strains and stresses at the Gauss Points.
		/// </summary>
		/// <param name="localDisplacements">The nodal displacements.</param>
		/// <returns>Two <see cref="IReadOnlyList{T}"/> with the strains and stresses of the Gauss points.</returns>
		public (IReadOnlyList<double[]> strains, IReadOnlyList<double[]> stresses)
			UpdateStrainStressesAtGaussPoints(double[] localDisplacements)
		{
			int numberOfGPs = QuadratureForStiffness.IntegrationPoints.Count;
			var strains = new double[numberOfGPs][];
			var stresses = new double[numberOfGPs][];
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < numberOfGPs; gp++)
			{
				IMatrixView constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				Matrix shapeGrandientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
				Matrix deformation = BuildDeformationMatrix(shapeGrandientsCartesian);

				strains[gp] = deformation.Multiply(localDisplacements);
				stresses[gp] = constitutive.Multiply(strains[gp]);
			}

			return (strains, stresses);
		}

		/// <summary>
		/// Assembles the deformation matrix of a solid element.
		/// The calculation are based on <see cref="https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch08.d/AFEM.Ch08.pdf"/>
		/// paragraph 8.4, equation 8.7.
		/// </summary>
		/// <param name="shapeGradientsCartesian">A <see cref="Matrix"/> containing the cartesian shape function gradients.</param>
		/// <returns>A <see cref="Matrix"/> containing the deformation matrix.</returns>
		private Matrix BuildDeformationMatrix(Matrix shapeGradientsCartesian)
		{
			var deformation = Matrix.CreateZero(6, 3 * Nodes.Count);
			for (int nodeIdx = 0; nodeIdx < Nodes.Count; nodeIdx++)
			{
				int col0 = 3 * nodeIdx;
				int col1 = 3 * nodeIdx + 1;
				int col2 = 3 * nodeIdx + 2;

				deformation[0, col0] = shapeGradientsCartesian[nodeIdx, 0];
				deformation[1, col1] = shapeGradientsCartesian[nodeIdx, 1];
				deformation[2, col2] = shapeGradientsCartesian[nodeIdx, 2];

				deformation[3, col0] = shapeGradientsCartesian[nodeIdx, 1];
				deformation[3, col1] = shapeGradientsCartesian[nodeIdx, 0];

				deformation[4, col1] = shapeGradientsCartesian[nodeIdx, 2];
				deformation[4, col2] = shapeGradientsCartesian[nodeIdx, 1];

				deformation[5, col0] = shapeGradientsCartesian[nodeIdx, 2];
				deformation[5, col2] = shapeGradientsCartesian[nodeIdx, 0];
			}

			return deformation;
		}

		/// <summary>
		/// The shape function matrix is 2-by-2n, where n = is the number of shape functions. Row 0 corresponds to dof X, while
		/// row 1 to dof Y, etc.
		/// </summary>
		private Matrix BuildShapeFunctionMatrix(double[] shapeFunctions)
		{
			var shapeFunctionMatrix = Matrix.CreateZero(3, 3 * shapeFunctions.Length);
			for (int i = 0; i < shapeFunctions.Length; i++)
			{
				shapeFunctionMatrix[0, 3 * i] = shapeFunctions[i];
				shapeFunctionMatrix[1, 2 * i + 1] = shapeFunctions[i];
				shapeFunctionMatrix[2, 3 * i + 2] = shapeFunctions[i];
			}
			return shapeFunctionMatrix;
		}

	}
}
