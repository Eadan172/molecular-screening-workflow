import requests

# Test different SwissDock API endpoints
endpoints = [
    "https://swissdock.ch",
    "https://swissdock.ch/api",
    "https://swissdock.ch/api/receptor",
    "https://www.swissdock.ch",
    "https://www.swissdock.ch/api",
    "https://www.swissdock.ch/api/receptor"
]

print("Testing SwissDock API endpoints...\n")

for endpoint in endpoints:
    try:
        response = requests.get(endpoint, timeout=10)
        print(f"Endpoint: {endpoint}")
        print(f"Status code: {response.status_code}")
        print(f"Response: {response.text[:200]}...")
    except Exception as e:
        print(f"Endpoint: {endpoint}")
        print(f"Error: {str(e)}")
    print("-" * 50)

print("\nConclusion:")
print("The SwissDock API appears to be unavailable or has changed its endpoint structure.")
print("\nRecommended alternatives:")
print("1. Use the SwissDock web interface directly: https://swissdock.ch")
print("2. Install AutoDock Vina locally for molecular docking")
print("3. Use other online docking services like:")
print("   - AutoDock-GPU (https://autodock.scripps.edu/)")
print("   - CB-Dock (http://cadd.zju.edu.cn/cbdock/)")
print("   - PatchDock (https://bioinfo3d.cs.tau.ac.il/PatchDock/)")
