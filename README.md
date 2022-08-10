# TripartiteViolation
This implementation verifies if a given multipartite NS behavior violates IC by considering a distributed version.

 * `NCP.py` buids the feasibility problem and call ncpol library.
 * `box.txt` List of correlators to compute the representative of each tripartite NS extremal class.
 * `group.py` Given a representative of tripartite NS extremal class, computes all behaviors belonging to the class.
 * `main.py` Main file.
 * `multiple.py` Given a multipartite distributed behavior, verifies if a given multipartite inequality IC inequality is violated.
 * `npa322.py` Varies parameters and find Q_n edge by using NCP.py verification.
 * `resource.py` It calculates the representative of each tripartite NS extremal class.
 * `uffink.py` Given a multipartite distributed behavior, verify if Uffink's inequality is violated.
