import h5py
import sys

def list_hdf5_contents(file_path):
    """
    List all datasets and groups in the HDF5 file.
    
    Parameters:
    file_path (str): Path to the HDF5 file
    """
    try:
        with h5py.File(file_path, 'r') as file:
            def print_structure(name, obj):
                if isinstance(obj, h5py.Dataset):
                    print(f"Dataset: {name}")
                    print(f"  Shape: {obj.shape}")
                    print(f"  Dtype: {obj.dtype}")
                    # Print first few values if dataset is small
                    try:
                        print(f"  Preview: {obj[:]}")
                    except Exception:
                        print("  (Unable to preview entire dataset)")
                elif isinstance(obj, h5py.Group):
                    print(f"Group: {name}")
            
            file.visititems(print_structure)
    
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except Exception as e:
        print(f"Error reading file: {e}")

def main():
    # Check if filename is provided as a command-line argument
    if len(sys.argv) < 2:
        print("Usage: python script.py <hdf5_filename>")
        sys.exit(1)
    
    filename = sys.argv[1]
    list_hdf5_contents(filename)

if __name__ == "__main__":
    main()
    