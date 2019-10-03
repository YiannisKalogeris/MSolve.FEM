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
	/// A two-dimensional truss element.
	/// </summary>
	public class Rod2D : IStructuralFiniteElement
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[2] { StructuralDof.TranslationX, StructuralDof.TranslationY };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private readonly double youngModulus;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		/// <summary>
		/// The truss density.
		/// </summary>
		public double Density { get; set; }

		/// <summary>
		/// The truss section.
		/// </summary>
		public double SectionArea { get; set; }


		/// <summary>
		/// Defines a <see cref="Rod2D"/> element.
		/// </summary>
		/// <param name="youngModulus">The truss material young modulus.</param>
		public Rod2D(double youngModulus)
		{
			this.youngModulus = youngModulus;
		}

		/// <summary>
		/// Defines a <see cref="Rod2D"/> element.
		/// </summary>
		/// <param name="youngModulus">The truss material young modulus.</param>
		/// <param name="dofEnumerator">An <see cref="IElementDofEnumerator"/>.</param>
		public Rod2D(double youngModulus, IElementDofEnumerator dofEnumerator)
			: this(youngModulus)
		{
			this.dofEnumerator = dofEnumerator;
		}

		/// <summary>
		/// Retrieves the type of the finite element used.
		/// </summary>
		public CellType CellType { get; } = CellType.Line;

		/// <summary>
		/// Defines the way that elemental degrees of freedom will be enumerated.
		/// For further info see <see cref="IElementDofEnumerator"/>.
		/// </summary>
		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		/// <summary>
		/// Calculates the transformation matrix of the <see cref="Rod2D"/> element.
		/// </summary>
		/// <param name="element">An element of type <see cref="Rod2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/>.</returns>
		public IMatrix TransformationMatrix(Element element)
		{
			double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			double L = Math.Sqrt(x2 + y2);
			double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
			double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;

			var transformation = Matrix.CreateZero(2, 4);
			transformation[0, 0] = c;
			transformation[0, 1] = s;
			transformation[1, 2] = c;
			transformation[1, 3] = s;
			return transformation;
		}

		/// <summary>
		/// Calculates the axial stresses throughout the element.
		/// </summary>
		/// <param name="element"><An element of type <see cref="Rod2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="local_d_Displacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>The value of the axial stress.</returns>
		public double CalculateAxialStress(Element element, double[] localDisplacements, double[] local_d_Displacements)
		{
			double[] globalStresses = CalculateStresses(element, localDisplacements, local_d_Displacements).Item2; // item1 = strains
			IMatrix transformation = TransformationMatrix(element);
			double[] localStresses = transformation.Multiply(globalStresses); // In local natural system there are 2 dofs
																			  // If Stress1 = localStresses[1] > 0 => tension. Else compression
			return localStresses[1];
		}

		#region IElementType Members

		/// <summary>
		/// The element ID.
		/// </summary>
		public int ID => 1;

		/// <summary>
		/// Dimensionality of the element.
		/// </summary>
		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="Rod2D"/>.</param>
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
		/// <param name="element">>An element of type <see cref="Rod2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="Rod2D"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element)
		{
			double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			double L = Math.Sqrt(x2 + y2);
			double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
			double c2 = c * c;
			double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
			double s2 = s * s;
			double cs = c * s;
			double E = this.youngModulus;
			double A = SectionArea;

			return dofEnumerator.GetTransformedMatrix(
				Matrix.CreateFromArray(new double[,]
				{
					{A*E*c2/L, A*E*cs/L, -A*E*c2/L, -A*E*cs/L },
					{A*E*cs/L, A*E*s2/L, -A*E*cs/L, -A*E*s2/L },
					{-A*E*c2/L, -A*E*cs/L, A*E*c2/L, A*E*cs/L },
					{-A*E*cs/L, -A*E*s2/L, A*E*cs/L, A*E*s2/L }
				}));
		}

		/// <summary>
		/// Calculates the mass matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="Rod2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="Rod2D"/>.</returns>
		public IMatrix MassMatrix(IElement element)
		{
			double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
			double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
			double L = Math.Sqrt(x2 + y2);

			double totalMassOver2 = Density * SectionArea * L / 2.0;

			int order = 4;
			var lumpedMass = Matrix.CreateZero(order, order);
			for (int i = 0; i < order; ++i) lumpedMass[i, i] = totalMassOver2;
			return lumpedMass;
		}

		/// <summary>
		/// Calculates the damping matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="Rod2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of an <see cref="Rod2D"/>.</returns>
		public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="Beam2DCorotationalAbstract"/>.</param>
		/// <param name="local_Displacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="local_d_Displacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="Tuple{T1,T2}"/> of the stresses and strains of the element.</returns>
		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] local_Displacements,
			double[] local_d_Displacements)
		{
			// WARNING: 1) No strains are computed 2) localdDisplacements are not used.
			double[] strains = null;
			double[] forces = CalculateForces(element, local_Displacements, local_d_Displacements);
			double[] stresses = Array.ConvertAll(forces, x => x / SectionArea);
			return new Tuple<double[], double[]>(strains, stresses);
		}

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="Rod2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
		{
			return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
		}

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="Rod2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			IMatrix stiffness = StiffnessMatrix(element);
			return stiffness.Multiply(localdDisplacements);
		}

		/// <summary>
		/// Calculates the forces applies to an <see cref="Rod2D"/> due to <see cref="MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="rod"/>.</param>
		/// <param name="loads">A list of <see cref="MassAccelerationLoad"/>. For more info see <seealso cref="MassAccelerationLoad"/>.</param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			var accelerations = new double[4];
			IMatrix massMatrix = MassMatrix(element);

			int index = 0;
			foreach (MassAccelerationLoad load in loads)
				foreach (IDofType[] nodalDOFTypes in dofs)
					foreach (IDofType dofType in nodalDOFTypes)
					{
						if (dofType == load.DOF) accelerations[index] += load.Amount;
						index++;
					}

			return massMatrix.Multiply(accelerations);
		}

		/// <summary>
		/// Save the current material state of the element.
		/// </summary>
		public void SaveMaterialState() { }

		#endregion

		#region IFiniteElement Members

		/// <summary>
		/// Boolean denoting if the material of the element has been modified.
		/// </summary>
		public bool MaterialModified => false;

		/// <summary>
		/// Resets any saved material states of the element to its initial state.
		/// </summary>
		public void ResetMaterialModified() { }

		#endregion

		#region IFiniteElement Members

		/// <summary>
		/// Clear the material state of the element.
		/// </summary>
		public void ClearMaterialState() { }

		/// <summary>
		/// Clear any saved material stresses of the element.
		/// </summary>
		public void ClearMaterialStresses() => throw new NotImplementedException();

		#endregion
	}
}
