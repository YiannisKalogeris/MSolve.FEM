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
	/// Enum that defines the possible directions of spring.
	/// </summary>
	public enum SpringDirections
	{
		X = 0, Y, Z, XY, YZ, XZ, XYZ
	}

	/// <summary>
	/// Three-dimensional spring element.
	/// </summary>
	public class SpringDamper3D : IStructuralFiniteElement
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private readonly double springCoefficient, dampingCoefficient;
		private readonly SpringDirections springDirections, dampingDirections;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		/// <summary>
		/// The element ID.
		/// </summary>
		public int ID => 999;

		/// <summary>
		/// Retrieves the type of the finite element used.
		/// </summary>
		public CellType CellType { get; } = CellType.Unknown;

		/// <summary>
		/// Dimensionality of the element.
		/// </summary>
		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

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
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="SpringDamper3D"/>.</param>
		/// <returns>A <see cref="IReadOnlyList{T}"/> that contains a <see cref="IReadOnlyList{T}"/> of <see cref="IDofType"/> with degrees of freedom for each elemental <see cref="Node"/>.</returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			if (element == null) return dofs;

			var d = new List<IReadOnlyList<IDofType>>();
			foreach (var node in element.Nodes)
			{
				var nodeDofs = new List<IDofType>();
				nodeDofs.AddRange(nodalDOFTypes);
				d.Add(nodeDofs);
			}
			return d;
		}

		/// <summary>
		/// Boolean denoting if the material of the element has been modified.
		/// </summary>
		public bool MaterialModified => false;

		/// <summary>
		/// Defines a <see cref="SpringDamper3D"/>.
		/// </summary>
		/// <param name="springCoefficient">The spring coefficient.</param>
		/// <param name="dampingCoefficient">The dampring coefficient.</param>
		/// <param name="springDirections">The direction of the spring.</param>
		/// <param name="dampingDirections">The direction of the damper.</param>
		public SpringDamper3D(double springCoefficient, double dampingCoefficient, SpringDirections springDirections,
			SpringDirections dampingDirections)
		{
			this.springCoefficient = springCoefficient;
			this.dampingCoefficient = dampingCoefficient;
			this.springDirections = springDirections;
			this.dampingDirections = dampingDirections;
		}

		/// <summary>
		/// Defines a <see cref="SpringDamper3D"/>.
		/// </summary>
		/// <param name="springCoefficient">The spring coefficient.</param>
		/// <param name="dampingCoefficient">The dampring coefficient.</param>
		/// <param name="springDirections">The direction of the spring.</param>
		/// <param name="dampingDirections">The direction of the damper.</param>
		/// <param name="dofEnumerator">An <see cref="IElementDofEnumerator"/>.</param>
		public SpringDamper3D(double springCoefficient, double dampingCoefficient, SpringDirections springDirections,
			SpringDirections dampingDirections, IElementDofEnumerator dofEnumerator)
			: this(springCoefficient, dampingCoefficient, springDirections, dampingDirections)
		{
			this.dofEnumerator = dofEnumerator;
		}

		/// <summary>
		/// Calculates the stiffness matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="SpringDamper3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="Beam3DCorotationalAbstract"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element)
		{
			double x = (springDirections == SpringDirections.X || springDirections == SpringDirections.XY || springDirections == SpringDirections.XZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
			double y = (springDirections == SpringDirections.Y || springDirections == SpringDirections.XY || springDirections == SpringDirections.YZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
			double z = (springDirections == SpringDirections.Z || springDirections == SpringDirections.XZ || springDirections == SpringDirections.YZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
			return Matrix.CreateFromArray(new double[,]
				{
					{x, 0, 0, -x, 0, 0},
					{0, y, 0, 0, -y, 0},
					{0, 0, z, 0, 0, -z},
					{-x, 0, 0, x, 0, 0},
					{0, -y, 0, 0, y, 0},
					{0, 0, -z, 0, 0, z}
				}
			);
		}

		/// <summary>
		/// Calculates the mass matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="SpringDamper3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="SpringDamper3D"/>.</returns>
		public IMatrix MassMatrix(IElement element) => Matrix.CreateZero(6, 6);

		/// <summary>
		/// Calculates the damping matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="SpringDamper3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of an <see cref="SpringDamper3D"/>.</returns>
		public IMatrix DampingMatrix(IElement element)
		{
			double x = (dampingDirections == SpringDirections.X || dampingDirections == SpringDirections.XY || dampingDirections == SpringDirections.XZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;
			double y = (dampingDirections == SpringDirections.Y || dampingDirections == SpringDirections.XY || dampingDirections == SpringDirections.YZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;
			double z = (dampingDirections == SpringDirections.Z || dampingDirections == SpringDirections.XZ || dampingDirections == SpringDirections.YZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;

			return Matrix.CreateFromArray(new double[,]
				{
				   { x, 0, 0, -x, 0, 0 },
				   { 0, y, 0, 0, -y, 0 },
				   { 0, 0, z, 0, 0, -z },
				   {-x, 0, 0, x, 0, 0 },
				   { 0,-y, 0, 0, y, 0 },
				   { 0, 0,-z, 0, 0, z }
				}
				);
		}

		/// <summary>
		/// Resets any saved material states of the element to its initial state.
		/// </summary>
		public void ResetMaterialModified() { }

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="SpringDamper3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="Tuple{T1,T2}"/> of the stresses and strains of the element.</returns>
		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements,
			double[] localdDisplacements)
			=> new Tuple<double[], double[]>(new double[6], new double[6]);

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="SpringDamper3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
			=> CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="SpringDamper3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			IMatrix stiffnessMatrix = StiffnessMatrix(element);
			return stiffnessMatrix.Multiply(localDisplacements);
		}

		/// <summary>
		/// Calculates the forces applies to an <see cref="SpringDamper3D"/> due to <see cref="MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="SpringDamper3D"/>.</param>
		/// <param name="loads">A list of <see cref="MassAccelerationLoad"/>. For more info see <seealso cref="MassAccelerationLoad"/>.</param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) => new double[6];

		/// <summary>
		/// Clear the material state of the element.
		/// </summary>
		public void ClearMaterialState() { }

		/// <summary>
		/// Save the current material state of the element.
		/// </summary>
		public void SaveMaterialState() { }

		/// <summary>
		/// Clear any saved material stresses of the element.
		/// </summary>
		public void ClearMaterialStresses() { }
	}
}
