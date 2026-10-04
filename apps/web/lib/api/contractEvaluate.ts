export function evaluateApiContract(
  expectedHash: string,
  runningHash: string,
): { status: "ok" | "mismatch"; scienceLocked: boolean } {
  const match = runningHash === expectedHash && /^[a-f0-9]{64}$/.test(runningHash);
  return {
    status: match ? "ok" : "mismatch",
    scienceLocked: !match,
  };
}
