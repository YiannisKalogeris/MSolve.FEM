using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.FreedomDegrees;
using MGroup.MSolve.Discretization.Interfaces;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Elements
{
	/// <summary>
	/// A beam of the euler-Bernoulli theory.
	/// </summary>
	public class EulerBeam2D : IStructuralFiniteElement
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.RotationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private readonly double youngModulus;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		/// <summary>
		/// The density of the beam.
		/// </summary>
		public double Density { get; set; }

		/// <summary>
		/// The section area of the beam.
		/// </summary>
		public double SectionArea { get; set; }

		/// <summary>
		/// The moment of inertia of the beam.
		/// </summary>
		public double MomentOfInertia { get; set; }

		/// <summary>
		/// The Raileigh alpha damping coefficient.
		/// </summary>
		public double RayleighAlpha { get; set; }

		/// <summary>
		/// The Raileigh beta damping coefficient.
		/// </summary>
		public double RayleighBeta { get; set; }

		/// <summary>
		/// Defines an <see cref="EulerBeam2D"/>.
		/// </summary>
		/// <param name="youngModulus">The young modulus of the beam.</param>
		public EulerBeam2D(double youngModulus)
		{
			this.youngModulus = youngModulus;
		}

		/// <summary>
		/// Defines an <see cref="EulerBeam2D"/>.
		/// </summary>
		/// <param name="youngModulus">The young modulus of the beam.</param>
		/// <param name="dofEnumerator">An <see cref="IElementDofEnumerator"/>.</param>
		public EulerBeam2D(double youngModulus, IElementDofEnumerator dofEnumerator)
			: this(youngModulus)
		{
			this.dofEnumerator = dofEnumerator;
		}

		/// <summary>
		/// Defines the way that elemental degrees of freedom will be enumerated.
		/// For further info see <see cref="IElementDofEnumerator"/>.
		/// </summary>
		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		#region IElementType Members

		/// <summary>
		/// The element ID.
		/// </summary>
		public int ID => 1;

		/// <summary>
		/// Retrieves the type of the finite element used.
		/// </summary>
		public CellType CellType { get; } = CellType.Line;

		/// <summary>
		/// Dimensionality of the element.
		/// </summary>
		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="EulerBeam2D"/>.</param>
		/// <returns>A <see cref="IReadOnlyList{T}"/> that contains a <see cref="IReadOnlyList{T}"/> of <see cref="IDofType"/> with degrees of freedom for each elemental <see cref="Node"/>.</returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofs;

		/// <summary>
		/// Retrieves the nodes of the element.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="EulerBeam2D"/>.</param>
		/// <returns>A list of the element nodes.</returns>
		public IList<Node> GetNodesForMatrixAssembly(Element element) => element.Nodes;

		/// <summary>
		/// Calculates the stiffness matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="EulerBeam2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="EulerBeam2D"/>.</returns>
		public virtual IMatrix StiffnessMatrix(IElement element)
		{
			double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			double L = Math.Sqrt(x2 + y2);
			double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
			double c2 = c * c;
			double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
			double s2 = s * s;
			double EL = this.youngModulus / L;
			double EAL = EL * SectionArea;
			double EIL = EL * MomentOfInertia;
			double EIL2 = EIL / L;
			double EIL3 = EIL2 / L;

			//TODO: optimize this
			int order = 6;
			var k = SymmetricMatrix.CreateFromPackedRowMajorArray(new double[]
			{
				c2*EAL+12*s2*EIL3, c*s*EAL-12*c*s*EIL3, -6*s*EIL2, -c2*EAL-12*s2*EIL3, -c*s*EAL+12*c*s*EIL3, -6*s*EIL2,
				s2*EAL+12*c2*EIL3, 6*c*EIL2, -s*c*EAL+12*c*s*EIL3, -s2*EAL-12*c2*EIL3, 6*c*EIL2,
				4*EIL, 6*s*EIL2, -6*c*EIL2, 2*EIL,
				c2*EAL+12*s2*EIL3, s*c*EAL-12*c*s*EIL3, 6*s*EIL2,
				s2*EAL+12*c2*EIL3, -6*c*EIL2,
				4*EIL
			}, order);

			return dofEnumerator.GetTransformedMatrix(k);
		}

		/// <summary>
		/// Calculates the mass matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="EulerBeam2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="EulerBeam2D"/>.</returns>
		public IMatrix MassMatrix(IElement element)
		{
			double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			double L = Math.Sqrt(x2 + y2);
			double L2 = L * L;
			double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
			double c2 = c * c;
			double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
			double s2 = s * s;
			double dAL420 = Density * SectionArea * L / 420;

			double totalMass = Density * SectionArea * L;
			double totalMassOfDiagonalTerms = 2 * dAL420 * (140 * c2 + 156 * s2) + 2 * dAL420 * (140 * s2 + 156 * c2);
			double scale = totalMass / totalMassOfDiagonalTerms;

			int order = 6;
			return SymmetricMatrix.CreateFromPackedRowMajorArray(new double[]
			{
				dAL420 *(140*c2+156*s2)*scale, 0, 0, 0, 0, 0,
				dAL420*(140*s2+156*c2)*scale, 0, 0, 0, 0,
				0, 0, 0, 0,
				dAL420*(140*c2+156*s2)*scale, 0, 0,
				dAL420*(140*s2+156*c2)*scale, 0,
				0
			}, order);
		}

		/// <summary>
		/// Calculates the damping matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="EulerBeam2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of an <see cref="EulerBeam2D"/>.</returns>
		public IMatrix DampingMatrix(IElement element)
		{
			var k = StiffnessMatrix(element);
			var m = MassMatrix(element);
			k.LinearCombinationIntoThis(RayleighBeta, m, RayleighAlpha);
			return k;
		}

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="EulerBeam2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="Tuple{T1,T2}"/> of the stresses and strains of the element.</returns>
		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements)
			=> throw new NotImplementedException();

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="EulerBeam2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
			=> CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="EulerBeam2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
			=> throw new NotImplementedException();


		/// <summary>
		/// Calculates the forces applies to an <see cref="EulerBeam2D"/> due to <see cref="MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="EulerBeam2D"/>.</param>
		/// <param name="loads">A list of <see cref="MassAccelerationLoad"/>. For more info see <seealso cref="MassAccelerationLoad"/>.</param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			var accelerations = new double[6];

			int index = 0;
			foreach (MassAccelerationLoad load in loads)
				foreach (IDofType[] nodalDOFTypes in dofs)
					foreach (IDofType dofType in nodalDOFTypes)
					{
						if (dofType == load.DOF) accelerations[index] += load.Amount;
						index++;
					}

			IMatrix massMatrix = MassMatrix(element);
			return massMatrix.Multiply(accelerations);
		}

		/// <summary>
		/// Save the current material state of the element.
		/// </summary>
		public void SaveMaterialState() => throw new NotImplementedException();

		#endregion

		#region IFiniteElement Members

		/// <summary>
		/// Boolean denoting if the material of the element has been modified.
		/// </summary>
		public bool MaterialModified => false;

		/// <summary>
		/// Resets any saved material states of the element to its initial state.
		/// </summary>
		public void ResetMaterialModified()
		{
			// Method intentionally left empty.
		}

		/// <summary>
		/// Clear the material state of the element.
		/// </summary>
		public void ClearMaterialState()
		{
			// Method intentionally left empty.
		}

		/// <summary>
		/// Clear any saved material stresses of the element.
		/// </summary>
		public void ClearMaterialStresses() => throw new NotImplementedException();
		#endregion
	}
}
