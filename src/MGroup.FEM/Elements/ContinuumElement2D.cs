﻿using System;
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

//TODO: Damping matrix calculation needs redesign for all of MSolve. For this class, see DampingMatrix().
//TODO: Materials also need redesign. Some properties are the same for all instances of a material class, some are the same for
//      all Gauss points of an element but differ across elements, some differ per Gauss point but are constant, and others
//      differ per Gauss point and per analysis iteration.
//TODO: Why does http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf Fig 31.9 take half the thickness when computing
//      the consistent mass matrix of Quad4 elements? For Tri3, the full thickness is used, as seen in Fig 31.7
//TODO: Simple Tri3 elements are more efficient than isoparamateric Tri3 elements. Many optimizations could be made. See
//      https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf for more.
//TODO: The shape of finite elements needs to be checked before the analysis: wrong node order, clockwise node order, the  
//      element's shape is too distorted, midpoints are too close to corners in quadratic elements, etc.
namespace MGroup.FEM.Elements
{
	/// <summary>
	/// Represents a continuum finite element for 2D problems. Specific elements (e.g. Quad4, Tri6, ...) can be created using
	/// the appropriate <see cref="IIsoparametricInterpolation2D"/>, <see cref="IQuadrature2D"/> etc. strategies. The thickness
	/// of this element is uniform, therefore it is necessary to use finer meshes to simulate domains with variable thickness.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class ContinuumElement2D : IStructuralFiniteElement, ICell<Node>
	{
		private readonly static IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY };
		private readonly IDofType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
		private DynamicMaterial dynamicProperties;
		private readonly IReadOnlyList<ElasticMaterial2D> materialsAtGaussPoints;

		public ContinuumElement2D(double thickness, IReadOnlyList<Node> nodes, IIsoparametricInterpolation2D interpolation,
			IQuadrature2D quadratureForStiffness, IQuadrature2D quadratureForConsistentMass,
			IGaussPointExtrapolation2D gaussPointExtrapolation,
			IReadOnlyList<ElasticMaterial2D> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
		{
			this.dynamicProperties = dynamicProperties;
			this.materialsAtGaussPoints = materialsAtGaussPoints;
			this.GaussPointExtrapolation = gaussPointExtrapolation;
			this.Nodes = nodes;
			this.Interpolation = interpolation;
			this.QuadratureForConsistentMass = quadratureForConsistentMass;
			this.QuadratureForStiffness = quadratureForStiffness;
			this.Thickness = thickness;

			dofTypes = new IDofType[nodes.Count][];
			for (int i = 0; i < nodes.Count; ++i)
			{
				dofTypes[i] = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY };
			}
		}

		public CellType CellType => Interpolation.CellType;
		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();
		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public int ID => throw new NotImplementedException(
			"Element type codes should be in a settings class. Even then it's a bad design choice");

		public IGaussPointExtrapolation2D GaussPointExtrapolation { get; }
		public IIsoparametricInterpolation2D Interpolation { get; }

		public bool MaterialModified
		{
			get
			{
				foreach (ElasticMaterial2D material in materialsAtGaussPoints)
					if (material.Modified) return true;
				return false;
			}
		}

		public IReadOnlyList<Node> Nodes { get; }
		public IQuadrature2D QuadratureForConsistentMass { get; }
		public IQuadrature2D QuadratureForStiffness { get; }
		public double Thickness { get; }

		public Matrix BuildConsistentMassMatrix()
		{
			int numDofs = 2 * Nodes.Count;
			var mass = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				//Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
				Matrix partial = shapeFunctionMatrix.MultiplyRight(shapeFunctionMatrix, true, false);
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				mass.AxpyIntoThis(partial, dA);
			}

			//WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
			mass.ScaleIntoThis(Thickness * dynamicProperties.Density);
			return mass;
		}

		public Matrix BuildLumpedMassMatrix()
		{
			int numDofs = 2 * Nodes.Count;
			var lumpedMass = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			// Contribution of each Gauss point to the element's area
			//TODO: Perhaps I could calculate the volume of the element without going through each Gauss point. Probably the 
			//      nodes are needed instead of the GPs. For linear elements I can find the area geometrically (as polygons).
			//TODO: this should have been cached when integrating other quantities (e.g. stiffness)
			double area = 0;
			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				area += jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
			}

			// Divide the total mass uniformly for each node
			double nodalMass = Thickness * area * dynamicProperties.Density / Nodes.Count;
			for (int i = 0; i < numDofs; ++i) lumpedMass[i, i] = nodalMass;

			return lumpedMass;
		}

		public Matrix BuildStiffnessMatrix()
		{
			int numDofs = 2 * Nodes.Count;
			var stiffness = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				// Calculate the necessary quantities for the integration
				IMatrixView constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				Matrix shapeGradientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
				Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);

				// Contribution of this gauss point to the element stiffness matrix
				Matrix partial = deformation.ThisTransposeTimesOtherTimesThis(constitutive);
				double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
				stiffness.AxpyIntoThis(partial, dA);
			}
			stiffness.ScaleIntoThis(Thickness);
			return stiffness;
		}

		//TODO: I think this method must be removed from IFiniteElement altogether. This procedure shoud be done for the global 
		//      mass matrix, once at the start of the dynamic analysis. The resulting vectors for each direction of the ground 
		//      motion should be stored. Then at each timestep they only need to be scaled and added to the rhs vector. The mass 
		//      matrix doesn't change, so there is not reason to recompute it at each time step.
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			int numDofs = 2 * Nodes.Count;
			var accelerations = new double[numDofs];
			IMatrix massMatrix = MassMatrix(element);

			foreach (MassAccelerationLoad load in loads)
			{
				int index = 0;
				foreach (IDofType[] nodalDOFTypes in dofTypes)
					foreach (IDofType dofType in nodalDOFTypes)
					{
						if (dofType == load.DOF) accelerations[index] += load.Amount;
						index++;
					}
			}

			return massMatrix.Multiply(accelerations);
		}

		public double CalculateArea()
		{
			//TODO: Linear elements can use the more efficient rules for volume of polygons. Therefore this method should be 
			//      delegated to the interpolation.
			//TODO: A different integration rule should be used for integrating constant functions. For linear elements there
			//      is only 1 Gauss point (most probably), therefore the computational cost could be the same as using the 
			//      polygonal formulas.
			double area = 0.0;
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				area += jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
			}
			return area;
		}

		public double[] CalculateForces(IElement element, double[] localTotalDisplacements, double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		//TODO: this method is probably not necessary at all
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
		{
			return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
		}

		//TODO: this method must be changed. It should calculates strains, stresses at GPs or nodes.
		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements,
			double[] localdDisplacements)
		{
			throw new NotImplementedException();
		}

		public void ClearMaterialState()
		{
			foreach (ElasticMaterial2D m in materialsAtGaussPoints) m.ClearState();
		}

		public void ClearMaterialStresses()
		{
			foreach (ElasticMaterial2D m in materialsAtGaussPoints) m.ClearStresses();
		}

		public IMatrix DampingMatrix(IElement element)
		{
			//TODO: Stiffness and mass matrices have already been computed probably. Reuse them.
			//TODO: Perhaps with Rayleigh damping, the global damping matrix should be created directly from global mass and stiffness matrices.
			Matrix damping = BuildStiffnessMatrix();
			damping.ScaleIntoThis(dynamicProperties.RayleighCoeffStiffness);
			damping.AxpyIntoThis(MassMatrix(element), dynamicProperties.RayleighCoeffMass);
			return damping;
		}

		/// <summary>
		/// Calculates the coordinates of the centroid of this element.
		/// </summary>
		public CartesianPoint FindCentroid()
			=> Interpolation.TransformNaturalToCartesian(Nodes, new NaturalPoint(0.0, 0.0));

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

		/// <summary>
		/// The returned structure is a list with as many entries as the number of nodes of this element. Each entry contains 
		/// a list with the dofs of the corresponding node. E.g. For node idx = 3, dof idx = 2 the IDof is result[3][2].
		/// </summary>
		/// <returns></returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetNodalDofs()
		{
			var allDofs = new IDofType[Nodes.Count][];
			for (int i = 0; i < Nodes.Count; ++i)
			{
				var nodalDofs = new IDofType[2];
				nodalDofs[0] = StructuralDof.TranslationX;
				nodalDofs[1] = StructuralDof.TranslationY;
				allDofs[i] = nodalDofs;
			}
			return allDofs;
		}

		//TODO: This is for the case when we also number constrained dofs globally.
		//// Perhaps this should be more minimalistic
		//// TODO: Either keep this or the GetNodlDofs() logic, but through a DofEnumerator and caching the dofs.
		//public DofTable<IDof> GetNodalDofsTable() 
		//{
		//    var dofTable = new DofTable<IDof>();
		//    int dofCounter = 0;
		//    foreach (var node in Nodes)
		//    {
		//        dofTable[node, DisplacementDof.X] = dofCounter++;
		//        dofTable[node, DisplacementDof.Y] = dofCounter++;
		//    }
		//    return dofTable;
		//}

		public IMatrix MassMatrix(IElement element)
		{
			return BuildConsistentMassMatrix();
			//return BuildLumpedMassMatrix();
		}

		public void ResetMaterialModified()
		{
			foreach (ElasticMaterial2D material in materialsAtGaussPoints) material.ResetModified();
		}

		public void SaveMaterialState()
		{
			foreach (ElasticMaterial2D m in materialsAtGaussPoints) m.SaveState();
		}

		public IMatrix StiffnessMatrix(IElement element) => DofEnumerator.GetTransformedMatrix(BuildStiffnessMatrix());

		/// <summary>
		/// Calculate strains (exx, eyy, 2exy) and stresses (sxx, syy, sxy) at integration points, store them in the materials 
		/// and return them (e.g. for postprocessing). The order of the tensors is the same as the order of the integration 
		/// points defined by <see cref="QuadratureForStiffness"/>.
		/// </summary>
		/// <param name="localDisplacements"></param>
		/// <returns></returns>
		public (IReadOnlyList<double[]> strains, IReadOnlyList<double[]> stresses) UpdateStrainsStressesAtGaussPoints(
			double[] localDisplacements)
		{
			int numGPs = QuadratureForStiffness.IntegrationPoints.Count;
			var strains = new double[numGPs][];
			var stresses = new double[numGPs][];
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < numGPs; ++gp)
			{
				IMatrixView constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				Matrix shapeGrandientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
				Matrix deformation = BuildDeformationMatrix(shapeGrandientsCartesian);

				strains[gp] = deformation.Multiply(localDisplacements);
				stresses[gp] = constitutive.Multiply(strains[gp]);
			}

			return (strains, stresses);
		}

		private Matrix BuildDeformationMatrix(Matrix shapeGradientsCartesian)
		{
			var deformation = Matrix.CreateZero(3, 2 * Nodes.Count);
			for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
			{
				int col0 = 2 * nodeIdx;
				int col1 = 2 * nodeIdx + 1;

				deformation[0, col0] = shapeGradientsCartesian[nodeIdx, 0];
				deformation[1, col1] = shapeGradientsCartesian[nodeIdx, 1];
				deformation[2, col0] = shapeGradientsCartesian[nodeIdx, 1];
				deformation[2, col1] = shapeGradientsCartesian[nodeIdx, 0];
			}
			return deformation;
		}

		/// <summary>
		/// The shape function matrix is 2-by-2n, where n = is the number of shape functions. Row 0 corresponds to dof X, while
		/// row 1 to dof Y.
		/// </summary>
		private Matrix BuildShapeFunctionMatrix(double[] shapeFunctions)
		{
			var shapeFunctionMatrix = Matrix.CreateZero(2, 2 * shapeFunctions.Length);
			for (int i = 0; i < shapeFunctions.Length; ++i)
			{
				shapeFunctionMatrix[0, 2 * i] = shapeFunctions[i];
				shapeFunctionMatrix[1, 2 * i + 1] = shapeFunctions[i];
			}
			return shapeFunctionMatrix;
		}
	}
}