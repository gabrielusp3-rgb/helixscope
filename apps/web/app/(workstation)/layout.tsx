import { WorkstationShell } from "@/components/shell/WorkstationShell";

export default function WorkstationLayout({ children }: { children: React.ReactNode }) {
  return <WorkstationShell>{children}</WorkstationShell>;
}
