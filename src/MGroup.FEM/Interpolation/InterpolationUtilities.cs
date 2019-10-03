using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.FEM.Interpolation
{
	/// <summary>
	/// A collection of methods that assist interpolations.
	/// Authors: Serafeim Bakalakos, Dimitris Tsapetis.
	/// </summary>
	public static class InterpolationUtilities
	{
		/// <summary>
		/// Transforms a natural point to cartesian.
		/// </summary>
		/// <param name="nodes">A <see cref="List{T}"/> of nodes.</param>
		/// <param name="shapeFunctionsAtNaturalPoint">The shape function at the natural point.</param>
		/// <returns>The <see cref="CartesianPoint"/> transformation of the <see cref="NaturalPoint"/>.</returns>
		public static CartesianPoint TransformPointNaturalToGlobalCartesian(IReadOnlyList<Node> nodes,
			Vector shapeFunctionsAtNaturalPoint)
		{
			int numFuncs = shapeFunctionsAtNaturalPoint.Length;
			if (nodes.Count != numFuncs) throw new ArgumentException(
				$"There are {numFuncs} evaluated shape functions stored, but {nodes.Count} were passed in.");
			double x = 0, y = 0;
			for (int i = 0; i < numFuncs; ++i)
			{
				Node node = nodes[i];
				double val = shapeFunctionsAtNaturalPoint[i];
				x += val * node.X;
				y += val * node.Y;
			}
			return new CartesianPoint(x, y);
		}

		/// <summary>
		/// Transforms a natural point to cartesian.
		/// </summary>
		/// <param name="nodes">A <see cref="List{T}"/> of nodes.</param>
		/// <param name="shapeFunctionsAtNaturalPoint">The shape function at the natural point.</param>
		/// <returns>The <see cref="CartesianPoint"/> transformation of the <see cref="NaturalPoint"/>.</returns>
		public static CartesianPoint TransformPointToGlobalCartesian(IReadOnlyList<Node> nodes,
			Vector shapeFunctionsAtNaturalPoint)
		{
			int numFuncs = shapeFunctionsAtNaturalPoint.Length;
			if (nodes.Count != numFuncs) throw new ArgumentException(
				$"There are {numFuncs} evaluated shape functions stored, but {nodes.Count} were passed in.");
			double x = 0, y = 0, z = 0;
			for (int i = 0; i < numFuncs; ++i)
			{
				Node node = nodes[i];
				double val = shapeFunctionsAtNaturalPoint[i];
				x += val * node.X;
				y += val * node.Y;
				z += val * node.Z;
			}

			return new CartesianPoint(x, y, z);
		}
	}
}
