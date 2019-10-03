namespace MGroup.FEM.Elements.SupportiveClasses
{
	/// <summary>
	/// An object containing the of a three dimensional integration point.
	/// </summary>
	public class ShapeFunctionNaturalDerivatives3D
	{
		#region Properties

		/// <summary>
		/// Shape function derivative per axis Heta.
		/// </summary>
		public double Eta { get; set; }

		/// <summary>
		/// Shape function derivative per axis Ksi.
		/// </summary>
		public double Xi { get; set; }

		/// <summary>
		/// Shape function derivative per axis Zeta.
		/// </summary>
		public double Zeta { get; set; }

		#endregion
	}
}
