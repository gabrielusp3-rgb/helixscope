"""OpenCL diagnosis: no fake devices, no random ICD injection."""

from __future__ import annotations

from modules import opencl_runtime


def test_hanging_gpu_icd_is_detected():
    assert opencl_runtime.ocl_icd_filenames_is_hanging(
        r"C:\Windows\System32\DriverStore\FileRepository\x\Intel_OpenCL_ICD64.dll"
    )
    assert not opencl_runtime.ocl_icd_filenames_is_hanging(
        r"C:\Program Files (x86)\Common Files\Intel\OpenCL\windows\intelocl64.dll"
    )


def test_subprocess_environ_strips_hanging_icd():
    env = opencl_runtime.subprocess_environ(
        {
            "OCL_ICD_FILENAMES": r"C:\Windows\System32\DriverStore\x\Intel_OpenCL_ICD64.dll",
            "PATH": "C:\\Windows",
        }
    )
    assert "OCL_ICD_FILENAMES" not in env or not opencl_runtime.ocl_icd_filenames_is_hanging(
        str(env.get("OCL_ICD_FILENAMES") or "")
    )
    assert env.get("HELIXSCOPE_OPENCL_ICD_STRIPPED") == "hanging_gpu_icd"


def test_device_selection_prefers_cpu_then_gpu():
    cpu = opencl_runtime.select_cas_offinder_device(
        {
            "n_platforms": 1,
            "platforms": [
                {
                    "devices": [
                        {"name": "Intel CPU", "type": "CPU"},
                        {"name": "Intel GPU", "type": "GPU"},
                    ]
                }
            ],
        }
    )
    assert cpu["device"] == "C"
    gpu = opencl_runtime.select_cas_offinder_device(
        {
            "n_platforms": 1,
            "platforms": [{"devices": [{"name": "Intel UHD", "type": "GPU"}]}],
        }
    )
    assert gpu["device"] == "G"
    none = opencl_runtime.select_cas_offinder_device({"n_platforms": 0, "platforms": []})
    assert none["device"] == "C"
    assert none["cpu_available"] is False
    assert opencl_runtime.cas_offinder_opencl_ready({"n_platforms": 0, "platforms": []}) is False
    assert opencl_runtime.cas_offinder_opencl_ready(
        {
            "n_platforms": 1,
            "platforms": [
                {
                    "name": "Intel(R) CPU Runtime",
                    "vendor": "Intel",
                    "version": "OpenCL 3.0",
                    "devices": [
                        {
                            "name": "Intel(R) Core",
                            "type": "CPU",
                            "vendor": "Intel",
                            "version": "OpenCL 3.0",
                        }
                    ],
                }
            ],
        }
    )


def test_diagnose_reports_blocker_or_ready():
    report = opencl_runtime.diagnose()
    assert "opencl_ready_for_cas_offinder" in report
    assert report["preferred_device"] in {"C", "G", "A"}
    assert report["gpu_required"] is False
    assert "official_cpu_runtime" in report
    if not report["opencl_ready_for_cas_offinder"]:
        assert report["blocker"]
        assert "OpenCL" in report["blocker"]
